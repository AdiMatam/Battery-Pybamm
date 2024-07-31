import numpy as np
import matplotlib.pyplot as plt

def generate_samples_uniform(mean, percent_diff, num_samples):
    """
    Generates random samples from a uniform distribution.
    
    Parameters:
    mean (float): The mean of the distribution.
    percent_diff (float): Percentage difference from the mean for bounds.
    num_samples (int): Number of samples to generate.
    
    Returns:
    numpy.ndarray: Array of samples from the uniform distribution.
    """
    lower_bound = mean - (percent_diff / 100) * mean
    upper_bound = mean + (percent_diff / 100) * mean
    samples = np.random.uniform(lower_bound, upper_bound, num_samples)
    return samples

def generate_samples_gaussian(mean, std_dev, num_samples):
    """
    Generates random samples from a Gaussian distribution.
    
    Parameters:
    mean (float): The mean of the distribution.
    std_dev (float): Standard deviation of the distribution.
    num_samples (int): Number of samples to generate.
    
    Returns:
    numpy.ndarray: Array of samples from the Gaussian distribution.
    """
    samples = np.random.normal(mean, std_dev, num_samples)
    return samples

# Parameters
mean = 0.385
percent_diff = 2
std_dev = 0.385 * 0.02 #0.005
num_samples = 1000  # Number of samples to generate

# Generate samples
uniform_samples = generate_samples_uniform(mean, percent_diff, num_samples)
gaussian_samples = generate_samples_gaussian(mean, std_dev, num_samples)

# Plotting histograms
plt.figure(figsize=(12, 6))

# Histogram for Uniform Distribution
plt.hist(uniform_samples, bins=30, alpha=0.5, color='blue', label='Uniform Distribution', density=True)

# Histogram for Gaussian Distribution
plt.hist(gaussian_samples, bins=30, alpha=0.5, color='red', label='Gaussian Distribution', density=True)

plt.xlabel('Sample Value')
plt.ylabel('Density')
plt.title('Histograms of Uniform and Gaussian Distribution Samples')
plt.legend()
plt.grid(True)
plt.show()
