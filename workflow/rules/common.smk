import polars as pl

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)

metadata = pl.read_csv(config['input_files']["metadata"], separator="\t", comment_prefix="#")

samples = metadata["Sample"].to_list()
read1 = metadata["R1"].to_list()
FIRST_SAMPLE = read1[0]
read2 = metadata["R2"].to_list()

READ_LOOKUP = dict(zip(samples, zip(read1, read2)))
