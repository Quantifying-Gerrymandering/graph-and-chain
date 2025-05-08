import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Your data
data = [50, 46, 45, 45, 45, 45, 45, 47, 45, 45, 47, 45, 45, 46, 45, 47, 45, 45, 46, 45, 47, 45, 45, 47, 45, 45, 46, 45, 47, 45, 47, 45, 45, 46, 45, 47, 45, 45, 46, 45, 43, 46, 47, 46, 45, 46, 45, 43, 46, 47, 46, 46, 46, 48, 46, 46, 47, 45, 46, 46, 45, 45, 46, 44, 46, 46, 44, 43, 47, 44, 46, 46, 47, 45]
data2 = [2, 6, 7, 7, 7, 7, 7, 5, 7, 7, 5, 7, 7, 6, 7, 5, 7, 7, 6, 7, 5, 7, 7, 5, 7, 7, 6, 7, 5, 7, 5, 7, 7, 6, 7, 5, 7, 7, 6, 7, 9, 6, 5, 6, 7, 6, 7, 9, 6, 5, 6, 6, 6, 4, 6, 6, 5, 7, 6, 6, 7, 7, 6, 8, 6, 6, 8, 9, 5, 8, 6, 6, 5, 7]

# Set up the figure
plt.figure(figsize=(10, 6))

# Create histogram
counts, bins, patches = plt.hist(data, bins=range(min(data), max(data)+2), 
                                  alpha=0.7, color='lightblue', edgecolor='darkblue',
                                  weights=np.ones(len(data)) / len(data) * 100)

# Calculate and plot kernel density estimation curve
x = np.linspace(min(data)-1, max(data)+1, 1000)

mu, std = stats.norm.fit(data)
pdf_values = stats.norm.pdf(x, mu, std)

bin_width = np.diff(bins)[0] 
y_normal = pdf_values * bin_width * 100

plt.plot(x, y_normal, 'darkblue', linewidth=2, label=f'Normal Fit: $\mu={mu:.2f}, \sigma={std:.2f}$')


# Add vertical line for observed value
observed_value = 46  # Example - replace with actual reference value
plt.axvline(x=observed_value, color='red', linewidth=2)

# Customize the plot
plt.xlabel('Democrat Seats Won')
plt.ylabel('Percent of Sampled Plans')
plt.title('Predicted Seats Won [FIX TITLE]')
plt.grid(axis='y', alpha=0.3)

# Set x-axis limits and ticks
plt.xlim(min(data)-2, max(data)+2)
plt.xticks(range(min(data)-2, max(data)+3))

# Add legend
plt.plot([], [], color='red', linewidth=2, label='CA2020')
plt.legend()

plt.tight_layout()
plt.show()