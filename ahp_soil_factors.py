import numpy as np

# Step 1: Create the pairwise comparison matrix
# Assume you have 4 criteria, and you've conducted your pairwise comparisons
# This is an example comparison matrix for illustration purposes
comparison_matrix = np.array([
    [1, 1/3,   1,   1/3,  1/3,    1/3,    1/3,    1/3,     1/3,    1/4,   1/5,    1/6,    1/5,   1/5,    1/3,    1/3,    1/3,     1/5,     1/4],
    [3, 1,     2,   1/2,  1/2,    1/2,    1/3,     2,       2,     1/2,   1/3,    1/3,    1/3,   1/3,    2,      2,      2,       1/2,     1/3],
    [1, 1/2,   1,   1/3,  1/3,    1/3,    1/3,    1/2,     1/3,    1/3,   1/4,    1/4,    1/3,   1/3,    1/2,    1/2,    1/2,     1/2,     1/3],
    [3, 2,     3,   1,    1/2,    1/2,    1/3,     1,       1/2,   1/2,   1/3,    1/3,    1/2,   1/2,    5,      2,      5,       1/2,     1/2],
    [3, 2,     3,   2,    1,       2,     1/2,     2,       2,     1/2,   1/2,    1/2,    1/2,   1/2,    2,      2,      2,       2,       1/2],
    [3, 2,     3,   2,    1/2,     1,     1/2,     2,       1/2,   1/2,   1/3,    1/3,    1/2,   1/2,     2,      2,      2,       2,       2],
    [3, 3,     3,   3,    2,       2,     1,       3,       5,     3,     1/3,    1/3,    1/3,   1/3,     3,      3,      3,       4,       3],
    [3, 1/2,   2,   1,    1/2,    1/2,    1/3,     1,       3,     1/3,   1/2,    1/2,    1/4,   1/4,     1,      1,      1,       2,       4],
    [3, 1/2,   3,   2,    1/2,    2,      1/5,     1/3,     1,     1/2,   1/4,    1/5,    1/5,   1/4,    2,      2,      2,       4,       2],
    [4, 2,     3,   2,    2,      2,      1/3,     3,       2,     1,     1/2,    1/4,    2,      2,     5,      5,      5,       4,       4],
    [5, 3,     4,   3,    2,      3,      3,       2,       4,     2,     1,      1/3,    2,     1/2,     5,      5,      5,       4,       4],
    [6, 3,     4,   3,    2,      3,      3,       2,       5,     4,     3,      1,      2,      2,     5,      5,      5,       5,       4],
    [5, 3,     3,   2,    2,      2,      3,       4,       5,     1/2,   1/2,    1/2,    1,     1/2,    6,      6,      6,       5,       2],
    [5, 3,     3,   2,    2,      2,      3,       4,       4,     1/2,   2,      1/2,    2,      1,     5,      5,      5,       4,       5],
    [3, 1/2,   2,   1,    1/2,    1/2,   1/3,     1/5,     1/2,    1/5,   1/5,    1/5,    1/6,   1/5,     1,      1/4,    4,       1/4,     1/4],
    [3, 1/2,   2,   1,    1/2,    1/2,   1/3,     1/2,     1/2,    1/5,   1/5,    1/5,    1/6,   1/5,    4,      1,      5,       1/2,     1/2],
    [3, 1/2,   2,   1,    1/2,    1/2,   1/3,     1/5,     1/2,    1/5,   1/5,    1/5,    1/6,   1/5,    1/4,    1/5,    1,       1/4,     1/4],
    [5,   2,    2,   2,    1/2,    1/2,   1/4,     1/2,     1/4,    1/4,   1/4,    1/5,    1/5,   1/4,    4,      2,      4,       1,       1/2],  
    [4,   3,    3,   2,    2,      1/2,   1/3,     1/4,     1/2,    1/4,   1/4,    1/4,    1/2,   1/5,    4,      2,      4,       2,       1]
])

# Step 2: Calculate the priority vector (weights) for the criteria
# Normalize the comparison matrix
column_sums = np.sum(comparison_matrix, axis=0)
normalized_matrix = comparison_matrix / column_sums

# Calculate the priority vector by taking the average across the rows
priority_vector = np.mean(normalized_matrix, axis=1)

# Step 3: Consistency Check
# Calculate the Consistency Index (CI)
lambda_max = np.mean(np.sum(comparison_matrix * priority_vector[:, np.newaxis], axis=0) / priority_vector)
CI = (lambda_max - len(priority_vector)) / (len(priority_vector) - 1)

# Calculate the Random Index (RI) for n=3 (you can find these values in AHP literature)
RI = 0.58  # Approximate value for n=3

# Calculate the Consistency Ratio (CR)
CR = CI / RI

print("Priority Vector (Weights):", priority_vector)
print("Consistency Ratio:", CR)

# Check the consistency
if CR < 0.1:
    print("The pairwise comparison matrix is consistent.")
else:
    print("The pairwise comparison matrix is not consistent. Please revise your comparisons.")
