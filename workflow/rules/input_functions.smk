# Input for Merging replicates 
def input_merge(wildcards):
    subset = samples.query("sample == @wildcards.sample")
    return {
        "r1": subset["R1"].tolist(),
        "r2": subset["R2"].tolist()
    }