import numpy as np
import random
import pandas as pd

# dy & dx distributions
p_dy = [0.8, 0.1, 0.05, 0.05]
print("p_dy distribution length & sum (should be 1):", len(p_dy), np.sum(p_dy))

p_dx = [0.1]*5 + [0.05]*10
print("p_dx distribution length & sum (should be 1):", len(p_dx), np.sum(p_dx))

num_cgi = 2752  # number of CGIs
cgi_dicts = []  # List of dictionaries to store each CGI row

for i in range(num_cgi):
    bR = random.randint(1500, 3000)  # CGI size limit
    print(f"\nLength of CGI {i+1} in bp:", bR)
    cgi = []
    l_cgi = 0

    while l_cgi < bR:
        # Sample dy
        r = random.uniform(0, 1)
        p_a, p_b = 0, p_dy[0]
        for j in range(len(p_dy)):
            if p_a < r <= p_b:
                dy = j * 2 + 2
                cgi.append(dy)
                l_cgi += dy
                break
            p_a += p_dy[j]
            if j < len(p_dy) - 1:
                p_b = p_a + p_dy[j+1]

        # Sample dx
        r = random.uniform(0, 1)
        p_a, p_b = 0, p_dx[0]
        for j in range(len(p_dx)):
            if p_a < r <= p_b:
                dx = j + 1
                cgi.append(dx)
                l_cgi += dx
                break
            p_a += p_dx[j]
            if j < len(p_dx) - 1:
                p_b = p_a + p_dx[j+1]

    # Add final dy
    r = random.uniform(0, 1)
    p_a, p_b = 0, p_dy[0]
    for j in range(len(p_dy)):
        if p_a < r <= p_b:
            dy = j * 2 + 2
            cgi.append(dy)
            l_cgi += dy
            break
        p_a += p_dy[j]
        if j < len(p_dy) - 1:
            p_b = p_a + p_dy[j+1]

    if len(cgi) % 2 == 0:
        print("Even-length CGI detected. Aborting this CGI. ................xxxxx.................")
    else:
        print(f"CGI {i+1} length in dx/dy values:", len(cgi))
        print(f"CGI {i+1} values:", cgi)

        # Convert to dictionary for DataFrame
        cgi_dict = {"Sequence": f"CGI_{i+1}"}
        for idx in range(0, len(cgi), 2):
            dy_col = f"dy{idx//2 + 1}"
            cgi_dict[dy_col] = cgi[idx]
            if idx + 1 < len(cgi):
                dx_col = f"dx{idx//2 + 1}"
                cgi_dict[dx_col] = cgi[idx + 1]
        cgi_dicts.append(cgi_dict)

# Create DataFrame
df = pd.DataFrame(cgi_dicts)

# Save to CSV
output_csv_path = r"F:\a_90percabove_excludedsequences\randomCGIs_matrix.csv"
df.to_csv(output_csv_path, index=False)
print(f"\nCGIs saved to: {output_csv_path}")
