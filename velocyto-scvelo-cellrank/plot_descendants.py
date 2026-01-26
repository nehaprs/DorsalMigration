import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/stage17_dm_fate_summary.csv")

plt.figure()
plt.bar(df["terminal_state"].values[:10], df["mean_fate_prob_stage1_interest"].values[:10])
plt.xticks(rotation=45, ha="right")
plt.ylabel("Mean fate probability")
plt.title("Stage1 cells of interest → terminal states")
plt.tight_layout()
plt.savefig("results/figures/stage17_dm_to_terminal_bar.png", dpi=200)


