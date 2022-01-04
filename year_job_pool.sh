#!/bin/bash

#set -e   # this doesn't work here for some reason
POOL_SIZE=3   # number of workers running in parallel

#######################################################################
#                            populate jobs                            #
#######################################################################

declare -a jobs

for (( i = 1983; i < 1990; i++ )); do
	jobs+=($i)
done

echo '################################################'
echo '    Launching jobs'
echo '################################################'

parallel() {
	local proc procs jobs cur
	jobs=("$@")         # input jobs array
	declare -a procs=() # processes array
	cur=0               # current job idx

	morework=true
	while $morework; do
		# if process array size < pool size, try forking a new proc
		if [[ "${#procs[@]}" -lt "$POOL_SIZE" ]]; then
			if [[ $cur -lt "${#jobs[@]}" ]]; then
				proc=${jobs[$cur]}
				echo "JOB ID = $cur; JOB = $proc."

				###############
				# do job here #
				###############

				tmpfile=./scripts/detect_ARs_polar_${proc}.py
				sed -e "98s/YEAR=[0-9]\{4\}/YEAR=$proc/" ./scripts/detect_ARs_polar.py > "$tmpfile"
				python "$tmpfile" &

				# add to current running processes
				procs+=("$!")
				# move to the next job
				((cur++))
			else
				morework=false
				continue
			fi
		fi

		for n in "${!procs[@]}"; do
			kill -0 "${procs[n]}" 2>/dev/null && continue
			# if process is not running anymore, remove from array
			unset procs[n]
		done
	done
	wait
}

parallel "${jobs[@]}"


