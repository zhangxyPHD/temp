#!/bin/bash 
#SBATCH -J We4.5-post
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=128
#SBATCH -o out1.txt    # Output file name
#SBATCH -e error1.txt  # Error file name
#SBATCH -t 70:00:00   # Time limit (hh:mm:ss)

# Parameter arrays
tsnap=0.01
Ldomain=8.0
tmax=20.0
Bos=(0.0)
epsilons=(0.01)
MAXlevels=(10)
Wes=(1 2 4 6 8 10 20 40 60 80 100)
Ohs=(0.001)
Js=(0 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2) #22

###################
# Parallel settings
###################
NPARA=32          # Maximum number of concurrent tasks
THREADS=4        # OpenMP threads per task

###################
# Compile
###################
cd basicmodel
if [ "$THREADS" -gt 1 ]; then
  qcc -w -Wall -O2 -disable-dimensions -fopenmp bounce.c -o bounce -lm
else
  qcc -w -Wall -O2 -disable-dimensions bounce.c -o bounce -lm
fi
qcc -w -Wall -O2 -disable-dimensions getFacet2D.c -o getFacet2D -lm
qcc -w -Wall -O2 -disable-dimensions getData2D-VP.c -o getData2D-VP -lm
cd ..

mkdir -p Results_Running

############################################
# Function: Check tasks and remove finished
############################################
# This function iterates over the 'tasks' array,
# removes any PID that has exited, and returns
# how many are still running.
check_and_clean_tasks() {
  local still_running=()
  for pid in "${tasks[@]}"; do
    if kill -0 "$pid" 2>/dev/null; then
      # PID is still alive
      still_running+=( "$pid" )
    fi
  done
  tasks=("${still_running[@]}")
  echo "${#tasks[@]}"  # Return how many remain
}

#################################
# Main loop over all parameters
#################################
declare -a tasks=()  # To store PIDs of launched jobs

for MAXlevel in "${MAXlevels[@]}"; do
  for We in "${Wes[@]}"; do
    for J in "${Js[@]}"; do
      for Oh in "${Ohs[@]}"; do
        for Bo in "${Bos[@]}"; do
          for epsilon in "${epsilons[@]}"; do

            FILENAME="Bo${Bo}-We${We}-J${J}-Oh${Oh}-MAXlevel${MAXlevel}-epsilon${epsilon}"
            CSVNAME="${FILENAME}.csv"

            # Skip if result already exists
            if [ -e "./Results_Running/$CSVNAME" ]; then
              echo "$FILENAME already exists."
              ##############################
              # Wait until concurrency < NPARA
              ##############################
              while true; do
                running_count=$(check_and_clean_tasks)
                if [ "$running_count" -lt "$NPARA" ]; then
                  break
                fi
                sleep 1
              done

              ###################################
              # Launch job in background
              ###################################
              (
                cd "$FILENAME" || exit 1
                if [ -e "video_$FILENAME.mp4" ]; then
                  echo "$FILENAME already exists."
                  continue
                fi
                export OMP_NUM_THREADS=1
                python3 getVideo.py --RMAX=3 --ZMAX=6 --tMAX=$tmax --tSNAP=$tsnap --CPUs=$THREADS >log_video 2>&1
                ffmpeg -framerate 30 -pattern_type glob -i 'Video/*.png' -vf scale=850:880 -c:v libx264 -r 30 -pix_fmt yuv420p video_$FILENAME.mp4 -y >log_video1 2>&1
              ) &

              # Record the PID of the background job
              tasks+=( "$!" )
            fi

          done
        done
      done
    done
  done
done

###############################
# Final wait for all tasks
###############################
# Check if any tasks still running, wait for them
while [ "${#tasks[@]}" -gt 0 ]; do
  # Wait for any job to finish
  wait -n 2>/dev/null || true
  # Clean up finished tasks from array
  check_and_clean_tasks >/dev/null
done

echo "All tasks have completed."
