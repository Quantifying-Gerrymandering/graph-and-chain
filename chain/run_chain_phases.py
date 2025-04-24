import geopandas as gpd
from gerrychain import Graph, Partition
import networkx as nx
import matplotlib.pyplot as plt
from functools import partial
import pandas as pd
import numpy as np
import math
import random
from gen_partition_random_starting_nodes import grow_districts
from gen_partition_spanning_tree import generate_spanning_tree_partition
from gerrychain.updaters import Tally, cut_edges
from gerrychain import MarkovChain
from gerrychain.constraints import contiguous, within_percent_of_ideal_population
from gerrychain.proposals import propose_random_flip
from gerrychain.optimization import SingleMetricOptimizer
from gerrychain.proposals import recom
from gerrychain.accept import always_accept
import json
from pathlib import Path

root = (Path(__file__).parent).parent # root directory of project
results = root / "chain" / "results"
final_partitions = root / "chain" / "chain-final-partitions"
final_scores = root / "chain" / "chain-final-scores"
final_stats = root / "chain" / "chain-final-stats"
e1_partitions = root / "chain" / "experiment1_initial_partitions"

# population metric: penalizes imbalances
def pop_with_comp(partition):
    avg_pop = np.average(list(partition["pop"].values()))
    return sum(pow(abs(pop - avg_pop), 5) for pop in list(partition["pop"].values()))

# compactness metric: number of cut edges
def comp(partition):
    return len(partition["cut_edges"])

# maximum pop deviation from the ideal
def pop_max_dev(partition):
    ideal_pop = sum(partition["pop"].values()) / len(partition["pop"])
    return max(abs(pop - ideal_pop) for pop in partition["pop"].values())

# linear piecewise beta function
def linear_piecewise_beta(hot_phases, cooldown_phases, cold_phases, phase_length):
    total_hot = hot_phases * phase_length
    total_cooldown = cooldown_phases * phase_length
    total_cold = cold_phases * phase_length
    def beta_function(step):
        if step < total_hot: 
            return 0
        elif step < total_hot + total_cooldown:
            discrete_phase = (step - total_hot) // phase_length
            return discrete_phase / cooldown_phases
        else:
            return 1
    return beta_function

# custom optimizer
class MyOptimizer(SingleMetricOptimizer):
    def _simulated_annealing_acceptance_function(self, beta_function, beta_magnitude):
        def acceptance_function(part):
            if part.parent is None: return True
            score_delta = self.score(part) - self.score(part.parent)
            beta = beta_function(part[self._step_indexer])
            if self._maximize: score_delta *= -1
            exponent = -beta * beta_magnitude * score_delta
            if exponent < -700: probability = 0.0
            elif exponent > 1: probability = 1.0
            else: probability = math.exp(exponent)
            return random.random() < probability
        return acceptance_function

# phase 1: optimize population (new implementation with additional steps to ensure population balance)
def run_phase1(chain_num, partition_file, phase_length, epsilon=0.05, initial_assignment=None):
    graph = Graph.from_json(partition_file)

    # Use provided initial assignment if given, otherwise use from file
    if initial_assignment is None:
        assignment_dict = {node: graph.nodes[node]["district_id"] for node in graph.nodes()}
    else:
        assignment_dict = initial_assignment
    
    initial_partition = Partition(
        graph,
        assignment=assignment_dict,
        updaters={
            "pop": Tally("CENS_Total", alias="pop"),
            "cut_edges": cut_edges,
        }
    )
    
    def calculate_max_deviation(partition):
        total_pop = sum(partition["pop"].values())
        num_districts = len(partition["pop"])
        ideal_pop = total_pop / num_districts
        max_dev = max(abs(pop - ideal_pop) / ideal_pop for pop in partition["pop"].values())
        return max_dev

    constraints = [contiguous]
    best_partition = None
    best_score = float("inf")
    worst_deviation = float("-inf")
    step_data = []
    
    optimizer = MyOptimizer(
        proposal=propose_random_flip,
        constraints=constraints,
        initial_state=initial_partition,
        optimization_metric=pop_with_comp,
        maximize=False
    )

    hot_phases, cooldown_phases, cold_phases = 10, 100, 40
    if chain_num == 0:
        beta_magnitude = 0.75
    elif chain_num == 1:
        beta_magnitude = 0.5
    else:
        beta_magnitude = 0.5

    # Initial number of steps
    total_steps = (hot_phases + cooldown_phases + cold_phases) * phase_length
    steps_completed = 0

    # run chain with different beta magnitudes for different initial partitions
    initial_chain = optimizer.simulated_annealing(
        total_steps,
        linear_piecewise_beta(hot_phases, cooldown_phases, cold_phases, phase_length),
        beta_magnitude=beta_magnitude,
        with_progress_bar=True
    )

    for i, part in enumerate(initial_chain):
        current_score = pop_with_comp(part)
        max_deviation = calculate_max_deviation(part)
        
        # Debug print for every 1000 steps
        if i % 1000 == 0:
            print(f"Step {i}: Max deviation = {max_deviation:.3%}")
        
        if current_score < best_score:
            best_score = current_score
            best_partition = part

        step_data.append({
            "iteration": steps_completed + i,
            "pop_score": current_score,
            "pop_max_dev": max_deviation,
            "comp": comp(part)
        })

    steps_completed += total_steps
    return best_partition, pd.DataFrame(step_data)

# custom implementation of within_percent_of_ideal_population
def custom_within_percent_of_ideal_population(epsilon, pop_key="pop"):
    def within_epsilon_of_ideal(partition):
        # calculate ideal population
        total_pop = sum(partition[pop_key].values())
        number_of_districts = len(partition[pop_key])
        ideal_pop = total_pop / number_of_districts
        
        # check population deviation
        for district_pop in partition[pop_key].values():
            # check if district population is within epsilon of ideal
            deviation = abs(district_pop - ideal_pop) / ideal_pop
            if deviation > epsilon:
                return False
        return True
    
    return within_epsilon_of_ideal

# phase 2: optimize compactness with a population constraint
def run_phase2(phase1_partition, partition_file, phase_length, epsilon=0.05):
    graph = Graph.from_json(partition_file)
    # convert assignment to dict
    if not isinstance(phase1_partition.assignment, dict): assignment_dict = dict(phase1_partition.assignment)
    else: assignment_dict = phase1_partition.assignment
    new_partition = Partition(
        graph,
        assignment=assignment_dict,
        updaters={
            "pop": Tally("CENS_Total", alias="pop"),
            "cut_edges": cut_edges,
        }
    )
    # calculate ideal population
    total_pop = sum(new_partition["pop"].values())
    num_districts = len(new_partition["pop"])
    ideal_pop = total_pop / num_districts
    print(f"Ideal population: {ideal_pop}")
    print(f"Number of districts: {num_districts}")

    hot_phases, cooldown_phases, cold_phases = 10, 100, 40
    total_steps = (hot_phases + cooldown_phases + cold_phases) * phase_length

    constraints = [
        contiguous,
        custom_within_percent_of_ideal_population(epsilon, pop_key="pop")
    ]

    optimizer = MyOptimizer(
        proposal=propose_random_flip,
        constraints=constraints,
        initial_state=new_partition,
        optimization_metric=comp,
        maximize=False
    )

    best_partition = None
    best_score = float("inf")
    step_data = []

    for i, part in enumerate(
        optimizer.simulated_annealing(
            total_steps,
            linear_piecewise_beta(hot_phases, cooldown_phases, cold_phases, phase_length),
            beta_magnitude=1.0,
            with_progress_bar=True
        )
    ):
        current_comp = comp(part)
        if current_comp < best_score:
            best_score = current_comp
            best_partition = part

        step_data.append({
            "iteration": i,
            "pop_score": pop_with_comp(part),
            "min_pop": min(part["pop"].values()),
            "max_pop": max(part["pop"].values()),
            "pop_max_dev": pop_max_dev(part),
            "comp_score": current_comp
        })

    return best_partition, pd.DataFrame(step_data)

def run_phase2_recom(phase1_partition, partition_file, phase_length, epsilon=0.05):
    graph = Graph.from_json(partition_file)
    
    # convert assignment to dict
    if not isinstance(phase1_partition.assignment, dict): assignment_dict = dict(phase1_partition.assignment)
    else: assignment_dict = phase1_partition.assignment

    initial_partition = Partition(
        graph,
        assignment=assignment_dict,
        updaters={
            "pop": Tally("CENS_Total", alias="pop"),
            "cut_edges": cut_edges,
        }
    )

    ideal_population = sum(initial_partition["pop"].values()) / len(initial_partition)

    proposal = partial(
        recom,
        pop_col="CENS_Total",
        pop_target=ideal_population,
        epsilon=epsilon,
        node_repeats=2
    )

    recom_chain = MarkovChain(
        proposal=proposal,
        constraints=[contiguous],
        accept=always_accept,
        initial_state=initial_partition,
        total_steps=phase_length
    )

    best_partition = None
    best_score = float("inf")
    step_data = []

    for i, part in enumerate(recom_chain):
        current_comp = comp(part)
        if current_comp < best_score:
            best_score = current_comp
            best_partition = part

        step_data.append({
            "iteration": i,
            "pop_score": pop_with_comp(part),
            "min_pop": min(part["pop"].values()),
            "max_pop": max(part["pop"].values()),
            "pop_max_dev": pop_max_dev(part),
            "comp_score": current_comp
        })

    print("Phase 2 done")
    print(f"Best compactness (cut edges): {best_score}")
    return best_partition, pd.DataFrame(step_data)

def calculate_file_id(chain_num, recom=False):
    # Calculates the file ID based on existing files.
    for i in range(1000):
        idstr = str(i).zfill(3)
        
        if chain_num == 0: path = f"spanning_tree_stats{idstr}"
        elif chain_num == 1: path = f"random_nodes_stats{idstr}"
        else: path = f"current_districting_stats{idstr}"

        if recom: path += "_recom.png"
        else: path += ".png"
        if not (final_stats/path).exists(): return str(i).zfill(3)
    return "999"

def run_chain(chain_num, partition_file, phase_length, epsilon=0.05, initial_assignment=None):
    # phase 1: optimize population and print stats
    phase1_best, phase1_stats = run_phase1(chain_num, partition_file, phase_length, epsilon, initial_assignment=initial_assignment)
    ideal_pop = sum(phase1_best["pop"].values()) / len(phase1_best["pop"])
    pop_array = np.array(list(phase1_best["pop"].values()))
    print("\nPhase 1 Stats:")
    print(f"   Min population: {pop_array.min()} ({abs((1 - pop_array.min()/ideal_pop)*100):.2f}% of ideal)")
    print(f"   Max population: {pop_array.max()} ({abs((1 - pop_array.max()/ideal_pop)*100):.2f}% of ideal)")

    # phase 2: optimize compactness with a population constraint
    epsilon = 0.05
    phase2_best, phase2_stats = run_phase2(phase1_best, partition_file, phase_length, epsilon)

    # phase 2: stats
    comp_array = comp(phase2_best)
    pop_array = np.array(list(phase2_best["pop"].values()))
    print("\nPhase 2 Stats:")
    print(f"   Compactness (cut edges): {comp_array}")
    print(f"   Min population: {pop_array.min()} ({abs((1 - pop_array.min()/ideal_pop)*100):.2f}% of ideal)")
    print(f"   Max population: {pop_array.max()} ({abs((1 - pop_array.max()/ideal_pop)*100):.2f}% of ideal)")
    combined_stats = pd.concat([phase1_stats, phase2_stats], ignore_index=True)

    # Save final scores and stats
    if chain_num == 0:
        scores_dict = {
            "phase1_pop_scores": phase1_stats['pop_score'].tolist(),
            "phase1_pop_max_dev": phase1_stats['pop_max_dev'].tolist(),
            "phase1_comp_scores": phase1_stats['comp'].tolist(),
            "phase2_pop_scores": phase2_stats['pop_score'].tolist(),
            "phase2_comp_scores": phase2_stats['comp_score'].tolist(),
            "phase2_pop_max_dev": phase2_stats['pop_max_dev'].tolist(),
            "all_iterations": combined_stats['iteration'].tolist()
        }
        with open(final_scores/"spanning_tree_scores.json", 'w') as f:
            json.dump(scores_dict, f)
            
    elif chain_num == 1:
        scores_dict = {
            "phase1_pop_scores": phase1_stats['pop_score'].tolist(),
            "phase1_pop_max_dev": phase1_stats['pop_max_dev'].tolist(),
            "phase1_comp": phase1_stats['comp'].tolist(),
            "phase2_pop_scores": phase2_stats['pop_score'].tolist(),
            "phase2_comp_scores": phase2_stats['comp_score'].tolist(),
            "phase2_pop_max_dev": phase2_stats['pop_max_dev'].tolist(),
            "all_iterations": combined_stats['iteration'].tolist()
        }

        with open(final_scores/"random_nodes_scores.json", 'w') as f:
            json.dump(scores_dict, f)
            
    else:
        scores_dict = {
            "phase1_pop_scores": phase1_stats['pop_score'].tolist(),
            "phase1_pop_max_dev": phase1_stats['pop_max_dev'].tolist(),
            "phase1_comp": phase1_stats['comp'].tolist(),
            "phase2_pop_scores": phase2_stats['pop_score'].tolist(),
            "phase2_comp_scores": phase2_stats['comp_score'].tolist(),
            "phase2_pop_max_dev": phase2_stats['pop_max_dev'].tolist(),
            "all_iterations": combined_stats['iteration'].tolist()
        }
        with open(final_scores/"current_districting_scores.json", 'w') as f:
            json.dump(scores_dict, f)

    # plot stats
    plt.figure(figsize=(12, 5))

    # phase 1: plot
    plt.subplot(1, 2, 1)
    plt.plot(phase1_stats["iteration"], phase1_stats["pop_score"], label="Population Metric")
    plt.plot(phase1_stats["iteration"], phase1_stats["comp"], label="Cut Edges")
    plt.title("Phase 1: Population Metric vs. Compactness Over Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Score")
    plt.legend()

    # phase 2: plot
    plt.subplot(1, 2, 2)
    plt.plot(phase2_stats["iteration"], phase2_stats["comp_score"], label="Compactness (Cut Edges)")
    plt.title("Phase 2: Compactness Over Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Cut Edges")
    plt.legend()

    plt.tight_layout()
    file_id = calculate_file_id(chain_num, recom=False)
    
    if chain_num == 0:
        with open(final_partitions/f"spanning_tree_final_partition{file_id}.json", 'w') as f:
            json.dump({str(k): v for k, v in phase2_best.assignment.items()}, f)
        plt.savefig(final_stats/f"spanning_tree_stats{file_id}.png")
    elif chain_num == 1:
        with open(final_partitions/f"random_nodes_final_partition{file_id}.json", 'w') as f:
            json.dump({str(k): v for k, v in phase2_best.assignment.items()}, f)
        plt.savefig(final_stats/f"random_nodes_stats{file_id}.png")
    else:
        with open(final_partitions/f"current_districting_final_partition{file_id}.json", 'w') as f:
            json.dump({str(k): v for k, v in phase2_best.assignment.items()}, f)
        plt.savefig(final_stats/f"current_districting_stats{file_id}.png")
    plt.close()
    return phase2_best, combined_stats

def run_chain_recom(chain_num, partition_file, phase_length, epsilon=0.05, initial_assignment=None):
    # phase 1: optimize population
    phase1_best, phase1_stats = run_phase1(chain_num, partition_file, phase_length, epsilon, initial_assignment=initial_assignment) # Pass through initial_assignment
    ideal_pop = sum(phase1_best["pop"].values()) / len(phase1_best["pop"])

    # phase 1: stats
    pop_array = np.array(list(phase1_best["pop"].values()))
    print("\nPhase 1 Stats:")
    print(f"   Min population: {pop_array.min()} ({abs((1 - pop_array.min()/ideal_pop)*100):.2f}% of ideal)")
    print(f"   Max population: {pop_array.max()} ({abs((1 - pop_array.max()/ideal_pop)*100):.2f}% of ideal)")

    # phase 2: optimize compactness with a population constraint
    epsilon = 0.05
    phase2_best, phase2_stats = run_phase2_recom(phase1_best, partition_file, phase_length, epsilon)

    # phase 2: stats
    comp_array = comp(phase2_best)
    pop_array = np.array(list(phase2_best["pop"].values()))
    print("\nPhase 2 Recom Stats:")
    print(f"   Compactness (cut edges): {comp_array}")
    print(f"   Min population: {pop_array.min()} ({abs((1 - pop_array.min()/ideal_pop)*100):.2f}% of ideal)")
    print(f"   Max population: {pop_array.max()} ({abs((1 - pop_array.max()/ideal_pop)*100):.2f}% of ideal)")

    # Combine phase1 and phase2 stats
    combined_stats = pd.concat([phase1_stats, phase2_stats], ignore_index=True)

    # Save final scores and stats
    if chain_num == 0:
        scores_dict = {
            "phase1_pop_scores": phase1_stats['pop_score'].tolist(),
            "phase1_pop_max_dev": phase1_stats['pop_max_dev'].tolist(),
            "phase1_comp_scores": phase1_stats['comp'].tolist(),
            "phase2_pop_scores": phase2_stats['pop_score'].tolist(),
            "phase2_comp_scores": phase2_stats['comp_score'].tolist(),
            "phase2_pop_max_dev": phase2_stats['pop_max_dev'].tolist(),
            "all_iterations": combined_stats['iteration'].tolist()
        }
        with open(final_scores/'spanning_tree_scores_recom.json', 'w') as f:
            json.dump(scores_dict, f)
            
    elif chain_num == 1:
        scores_dict = {
            "phase1_pop_scores": phase1_stats['pop_score'].tolist(),
            "phase1_pop_max_dev": phase1_stats['pop_max_dev'].tolist(),
            "phase1_comp": phase1_stats['comp'].tolist(),
            "phase2_pop_scores": phase2_stats['pop_score'].tolist(),
            "phase2_comp_scores": phase2_stats['comp_score'].tolist(),
            "phase2_pop_max_dev": phase2_stats['pop_max_dev'].tolist(),
            "all_iterations": combined_stats['iteration'].tolist()
        }
        with open(final_scores/'random_nodes_scores_recom.json', 'w') as f:
            json.dump(scores_dict, f)
            
    else:
        scores_dict = {
            "phase1_pop_scores": phase1_stats['pop_score'].tolist(),
            "phase1_pop_max_dev": phase1_stats['pop_max_dev'].tolist(),
            "phase1_comp": phase1_stats['comp'].tolist(),
            "phase2_pop_scores": phase2_stats['pop_score'].tolist(),
            "phase2_comp_scores": phase2_stats['comp_score'].tolist(),
            "phase2_pop_max_dev": phase2_stats['pop_max_dev'].tolist(),
            "all_iterations": combined_stats['iteration'].tolist()
        }
        with open(final_scores/'/current_districting_scores_recom.json', 'w') as f:
            json.dump(scores_dict, f)

    # plot stats
    plt.figure(figsize=(12, 5))

    # phase 1: plot
    plt.subplot(1, 2, 1)
    plt.plot(phase1_stats["iteration"], phase1_stats["pop_score"], label="Population Metric")
    plt.plot(phase1_stats["iteration"], phase1_stats["comp"], label="Cut Edges")
    plt.title("Phase 1: Population Metric vs. Compactness Over Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Score")
    plt.legend()

    # phase 2: plot
    plt.subplot(1, 2, 2)
    plt.plot(phase2_stats["iteration"], phase2_stats["comp_score"], label="Compactness (Cut Edges)")
    plt.title("Phase 2 Recom: Compactness Over Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Cut Edges")
    plt.legend()

    plt.tight_layout()
    file_id = calculate_file_id(chain_num, True)
    if chain_num == 0:
        with open(final_partitions/f'spanning_tree_final_partition{file_id}_recom.json', 'w') as f:
            json.dump({str(k): v for k, v in phase2_best.assignment.items()}, f)
        plt.savefig(final_stats/f"spanning_tree_stats{file_id}_recom.png")
    elif chain_num == 1:
        with open(final_partitions/f'random_nodes_final_partition{file_id}_recom.json', 'w') as f:
            json.dump({str(k): v for k, v in phase2_best.assignment.items()}, f)
        plt.savefig(final_stats/f"random_nodes_stats{file_id}_recom.png")
    else:
        with open(final_partitions/f'current_districting_final_partition{file_id}_recom.json', 'w') as f:
            json.dump({str(k): v for k, v in phase2_best.assignment.items()}, f)
        plt.savefig(final_stats/f"current_districting_stats{file_id}_recom.png")
    plt.close()
    return phase2_best, combined_stats

def verify_partition(graph, partition):
    """Verify that each district in the partition is contiguous."""
    districts = set(partition.values())
    print(f"\nChecking {len(districts)} districts for contiguity...")
    
    non_contiguous = []
    for district_id in districts:
        district_nodes = [node for node, dist in partition.items() if dist == district_id]
        subgraph = graph.subgraph(district_nodes)
        
        # Check if connected
        if not nx.is_connected(subgraph):
            non_contiguous.append(district_id)
            print(f"\nDistrict {district_id} is not contiguous!")
            components = list(nx.connected_components(subgraph))
            print(f"  Components: {len(components)}")
            print(f"  Sizes: {[len(c) for c in components]}")
            
            # Sort components by size (largest to smallest)
            components = sorted(components, key=len, reverse=True)
            main_component = components[0]
            
            # For each smaller component, find nearest node in main component
            for small_comp in components[1:]:
                print("\nDisconnected component GEOIDs:")
                for node in small_comp:
                    print(f"  {gdf.iloc[node]['GEOID20']}")
                
                # Find closest node in main component
                min_distance = float('inf')
                closest_pair = None
                for node1 in small_comp:
                    point1 = gdf.iloc[node1].geometry.centroid
                    for node2 in main_component:
                        point2 = gdf.iloc[node2].geometry.centroid
                        dist = point1.distance(point2)
                        if dist < min_distance:
                            min_distance = dist
                            closest_pair = (node1, node2)
                
                if closest_pair:
                    print("\nClosest connection:")
                    print(f"  Disconnected: {gdf.iloc[closest_pair[0]]['GEOID20']}")
                    print(f"  Main component: {gdf.iloc[closest_pair[1]]['GEOID20']}")
                    print(f"  Distance: {min_distance}")
    
    if non_contiguous:
        print(f"\nFound {len(non_contiguous)} non-contiguous districts: {non_contiguous}")
        return False
    else:
        print("All districts are contiguous!")
        return True

if __name__ == "__main__":
    partition_file = root/"graph"/"dual_graph.json"
    shapefile = root/"data"/"cacbg20"/"cacbg20.shp"
    
    hot, cooldown, cold = 10, 100, 40
    phase_len = 1000
    epsilon = 0.05
    graph = Graph.from_json(partition_file)
    gdf = gpd.read_file(shapefile)

    for i in range(0, 10):
        random_nodes_partition = grow_districts(graph, 52, gdf)
        run_chain(chain_num=1,
                  partition_file=partition_file, 
                  phase_length=phase_len, 
                  epsilon=epsilon, 
                  initial_assignment=random_nodes_partition,
                  file_id=i)