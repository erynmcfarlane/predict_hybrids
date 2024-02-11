#!/bin/bash
#SBATCH --job-name=dfuse-uneven
#SBATCH --nodes=1
#SBATCH --time=0-8:00:00
#SBATCH --account=modelscape
#SBATCH --mem-per-cpu=500G

module load arcc/1.0
module load gcc/12.2.0
module load gsl/2.7.1
module load r

### all new dfuse simulations to look at how assymetry affects consisency in SNPs that are or aren't under selection ####
### all of these are 3 deme cases, rather than 11 deme cases
### for the revision, Zach wants all of these in all 12 scenarios, not just two, so, here we go.


### 100:1

./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./hundredone/dmi_m0.01_neutral
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./hundredone/dmi_m0.01_c0.2
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./hundredone/dmi_m0.01_c0.9

./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./hundredone/dmi_m0.2_neutral
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./hundredone/dmi_m0.2_c0.2
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./hundredone/dmi_m0.2_c0.9

###path
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./hundredone/path_m0.01_neutral
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./hundredone/path_m0.01_c0.2
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./hundredone/path_m0.01_c0.9

./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./hundredone/path_m0.2_neutral
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./hundredone/path_m0.2_c0.2
./dfuse -d 3_deme_uneven_100_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./hundredone/path_m0.2_c0.9

### 10:1
### DMI

./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./tenone/dmi_m0.01_neutral
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./tenone/dmi_m0.01_c0.2
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./tenone/dmi_m0.01_c0.9

./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./tenone/dmi_m0.2_neutral
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./tenone/dmi_m0.2_c0.2
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./tenone/dmi_m0.2_c0.9

###path
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./tenone/path_m0.01_neutral
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./tenone/path_m0.01_c0.2
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./tenone/path_m0.01_c0.9

./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./tenone/path_m0.2_neutral
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./tenone/path_m0.2_c0.2
./dfuse -d 3_deme_uneven_10_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./tenone/path_m0.2_c0.9



./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./fiveone/dmi_m0.01_neutral
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./fiveone/dmi_m0.01_c0.2
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./fiveone/dmi_m0.01_c0.9

./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./fiveone/dmi_m0.2_neutral
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./fiveone/dmi_m0.2_c0.2
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./fiveone/dmi_m0.2_c0.9

###path
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./fiveone/path_m0.01_neutral
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./fiveone/path_m0.01_c0.2
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./fiveone/path_m0.01_c0.9

./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./fiveone/path_m0.2_neutral
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./fiveone/path_m0.2_c0.2
./dfuse -d 3_deme_uneven_5_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./fiveone/path_m0.2_c0.9


./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./twoone/dmi_m0.01_neutral
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./twoone/dmi_m0.01_c0.2
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./twoone/dmi_m0.01_c0.9

./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./twoone/dmi_m0.2_neutral
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./twoone/dmi_m0.2_c0.2
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_1.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./twoone/dmi_m0.2_c0.9

###path
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.01 -c 0 -o ./twoone/path_m0.01_neutral
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.2 -o ./twoone/path_m0.01_c0.2
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.01 -c 0.9 -o ./twoone/path_m0.01_c0.9

./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1_neutral.txt -r 20 -g 100 -m 0.2 -c 0 -o ./twoone/path_m0.2_neutral
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.2 -o ./twoone/path_m0.2_c0.2
./dfuse -d 3_deme_uneven_2_1.txt -e /gscratch/buerkle/data/incompatible/runs/inputfiles/epifile_2c.txt -s /gscratch/buerkle/data/incompatible/runs/inputfiles/selfile_1.txt -r 20 -g 100 -m 0.2 -c 0.9 -o ./twoone/path_m0.2_c0.9
