#parse samples
samples = pd.read_csv("config/samples.tsv", dtype=str, sep="\t").set_index(["sample"], drop=False)
SAMPLES = samples['sample'].to_list()
