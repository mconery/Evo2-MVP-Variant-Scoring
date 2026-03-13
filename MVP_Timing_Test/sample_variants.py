import pandas as pd
import os

input_file = "/vast/projects/anuragv/cohort/mconery/mvp_variant_test/MVP_matched_variants.csv"
output_dir = "/vast/projects/anuragv/cohort/mconery/mvp_timing_test"
output_file = os.path.join(output_dir, "timing_sampled_variants.csv")

df = pd.read_csv(input_file)

low_pip = df[df['Set'] == 'Low_PIP'].sample(50, random_state=42)
high_pip = df[df['Set'] == 'High_PIP'].sample(50, random_state=42)

sampled = pd.concat([low_pip, high_pip]).reset_index(drop=True)

os.makedirs(output_dir, exist_ok=True)
sampled.to_csv(output_file, index=False)
print(f"Saved {len(sampled)} variants to {output_file}")
