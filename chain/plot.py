import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Your data
data = [50, 46, 45, 45, 45, 45, 45, 47, 45, 45, 47, 45, 45, 46, 45, 47, 45, 45, 46, 45, 
        47, 45, 45, 47, 45, 45, 46, 45, 47, 45, 47, 45, 45, 46, 45, 47, 45, 45, 46, 45, 
        43, 46, 47, 46, 45, 46, 45, 43, 46, 47, 46]

# Set up the figure
plt.figure(figsize=(10, 6))

# Create histogram
counts, bins, patches = plt.hist(data, bins=range(min(data), max(data)+2), 
                                  alpha=0.7, color='lightblue', edgecolor='darkblue',
                                  weights=np.ones(len(data)) / len(data) * 100)

# Calculate and plot kernel density estimation curve
x = np.linspace(min(data)-1, max(data)+1, 1000)
kde = stats.gaussian_kde(data, bw_method=0.5)
y = kde(x) * np.diff(bins)[0] * 100  # Scale to match histogram percentage
plt.plot(x, y, 'darkblue', linewidth=2)

# Add vertical line for observed value
observed_value = 46  # Example - replace with actual reference value
plt.axvline(x=observed_value, color='red', linewidth=2)

# Customize the plot
plt.xlabel('Democrat Seats Won')
plt.ylabel('Percent of Sampled Plans')
plt.title('Duke University NC 2016 USH Vote Outcomes')
plt.grid(axis='y', alpha=0.3)

# Set x-axis limits and ticks
plt.xlim(40, 52)
plt.xticks(range(40, 53))

# Add legend
plt.plot([], [], color='red', linewidth=2, label='NC2016')
plt.legend()

plt.tight_layout()
plt.show()