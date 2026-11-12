import numpy as np
import random
import pandas as pd
import os


# RANDOM SEED


random.seed(42)
np.random.seed(42)


# PROBABILITY DISTRIBUTIONS


# dy values: 2, 4, 6, 8
p_dy = [0.8, 0.1, 0.05, 0.05]
dy_values = [2, 4, 6, 8]

# dx values: 1 to 15
p_dx = [0.1] * 5 + [0.05] * 10
dx_values = list(range(1, 16))

print(
    "p_dy distribution length and sum:",
    len(p_dy),
    sum(p_dy)
)

print(
    "p_dx distribution length and sum:",
    len(p_dx),
    sum(p_dx)
)


# GENERATION SETTINGS


num_cgi = 2752
num_sets = 100

MIN_CGI_LENGTH = 1500
MAX_CGI_LENGTH = 3000


# OUTPUT DIRECTORY


output_dir = (
    "path_to_the_output_folder"
)

os.makedirs(output_dir, exist_ok=True)


# FUNCTION: GENERATE ONE RANDOM CGI


def generate_random_cgi():

    # Random target length
    bR = random.randint(
        MIN_CGI_LENGTH,
        MAX_CGI_LENGTH
    )

    cgi = []
    l_cgi = 0


    # GENERATE ALTERNATING dy AND dx VALUES


    while l_cgi < bR:

        # Sample dy
        dy = random.choices(
            dy_values,
            weights=p_dy,
            k=1
        )[0]

        cgi.append(dy)
        l_cgi += dy

        # Sample dx
        dx = random.choices(
            dx_values,
            weights=p_dx,
            k=1
        )[0]

        cgi.append(dx)
        l_cgi += dx


    # ADD FINAL dy


    dy = random.choices(
        dy_values,
        weights=p_dy,
        k=1
    )[0]

    cgi.append(dy)
    l_cgi += dy


    # CHECK ALTERNATING STRUCTURE


    if len(cgi) % 2 == 0:

        raise ValueError(
            "Even-length CGI detected."
        )

    return cgi, l_cgi, bR



# GENERATE 100 INDEPENDENT RANDOM CGI SETS


for file_num in range(1, num_sets + 1):

    
    print(f"GENERATING RANDOM CGI SET {file_num}/100")
    

    cgi_dicts = []


    # GENERATE 2752 CGIs FOR THIS SET


    for i in range(num_cgi):

        cgi, actual_length, target_length = (
            generate_random_cgi()
        )


        # CONVERT CGI INTO DICTIONARY


        cgi_dict = {
            "Sequence": f"CGI_{i + 1}"
        }

        for idx in range(0, len(cgi), 2):

            dy_col = f"dy{idx // 2 + 1}"

            cgi_dict[dy_col] = cgi[idx]

            if idx + 1 < len(cgi):

                dx_col = f"dx{idx // 2 + 1}"

                cgi_dict[dx_col] = cgi[idx + 1]

        cgi_dicts.append(cgi_dict)


    # CREATE DATAFRAME


    df = pd.DataFrame(cgi_dicts)


    # SAVE THIS INDEPENDENT SET


    output_csv_path = os.path.join(
        output_dir,
        f"randomCGIs_matrix_{file_num}.csv"
    )

    df.to_csv(
        output_csv_path,
        index=False
    )

    print(
        f"Random CGI set {file_num}/100 saved."
    )

    print(
        "Number of CGIs:",
        len(df)
    )

    print(
        "Output:",
        output_csv_path
    )


# FINAL STATUS



print("RANDOM CGI GENERATION COMPLETED")


print("Total independent sets:", num_sets)
print("CGIs per set:", num_cgi)
print("Total generated CGIs:", num_sets * num_cgi)

print("\nOutput directory:")
print(output_dir)

print("\nDONE.")
