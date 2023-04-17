#!/bin/bash

########################################################### submitting params
partition="v6_384"
nodes=1
ntasks_per_node=10
cpus_per_task=1
mem_per_cpu=16G


########################################################### program params
# executable object
exe="/public1/home/sc81233/xuchengwang/pqmc-framework/build/pqmc_framework"

# program options
# edit the config file
config_file="/public1/home/sc81233/xuchengwang/pqmc-framework/benchmark/config.toml"

# model params
t=1.0
u=4.0
sed -i "s/ t = .*/ t = "$t"/" $config_file
sed -i "s/ u = .*/ u = "$u"/" $config_file

l=4
sed -i "s/ nl = .*/ nl = "$l"/" $config_file

# momentum for measurments of observables
momentum="MPoint"
momentum_list="KstarsAll"
sed -i "s/momentum = .*/momentum = \""$momentum"\"/" $config_file
sed -i "s/momentum_list = .*/momentum_list = \""$momentum_list"\"/" $config_file

# monte carlo params
# inverse temperature, imaginary-time grids and pace of stabilization
theta=5.0
dt=0.05
nt=200
beta=1.25
ntm=50
stabilization_pace=10
sed -i "s/ theta = .*/ theta = "$theta"/" $config_file
sed -i "s/ beta = .*/ beta = "$beta"/" $config_file
sed -i "s/ dt = .*/ dt = "$dt"/" $config_file
sed -i "s/ nt = .*/ nt = "$nt"/" $config_file
sed -i "s/ ntm = .*/ ntm = "$ntm"/" $config_file
sed -i "s/stabilization_pace = .*/stabilization_pace = "$stabilization_pace"/" $config_file

sweeps_warmup=1000
bin_num=2000
bin_capacity=50
sweeps_between_bins=50
sed -i "s/sweeps_warmup = .*/sweeps_warmup = "$sweeps_warmup"/" $config_file
sed -i "s/bin_num = .*/bin_num = "$bin_num"/" $config_file
sed -i "s/bin_capacity = .*/bin_capacity = "$bin_capacity"/" $config_file
sed -i "s/sweeps_between_bins = .*/sweeps_between_bins = "$sweeps_between_bins"/" $config_file

sed -i "s/observables = .*/observables = [ \"DoubleOccupation\", \"KineticEnergy\" ]/" $config_file

# set up name of the output folder 
folder_name="l"$l"u"$u
output_folder="/public1/home/sc81233/xuchengwang/pqmc-framework/benchmark/"$folder_name
ising_fields_file=$output_folder"/ising.fields"

# set up jobname and log output name
job_name=$folder_name
log_file=$output_folder"/log.out"
err_file=$output_folder"/err.out"

# create output folder if not exist
if [ ! -d $output_folder ]; then
    mkdir -p $output_folder
fi
# delete previous log file if exist
if [ -f $log_file ]; then
    rm $log_file && touch $log_file
fi


########################################################### submit mission
sbatch \
--job-name=$job_name --output=$log_file --error=$err_file --partition=$partition \
--nodes=$nodes --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task --mem-per-cpu=$mem_per_cpu \
--export=exe=$exe,config_file=$config_file,output_folder=$output_folder,ising_fields_file=$ising_fields_file \
./run.sh


exit 0