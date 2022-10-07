root_dir=/usr2/people/yolandahuang/pacbio/YHPB005
in_dir=$root_dir/fastq/split_files
out_dir=out_2207
python3 src/run_steps.py old_gibson_cfg.json $in_dir $out_dir 6
