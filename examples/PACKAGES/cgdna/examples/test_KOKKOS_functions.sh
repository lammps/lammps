###################################################################################################
# "Double Maths" Machinary
###################################################################################################

find_double_maths_files() {
	mapfile -t DOUBLE_MATHS_FILES < <(
		cd "$LMPDIR/src/KOKKOS" &&
			find . -type f -exec grep -l 'sqrtf(\|expf(' {} + | sort
	)
}

ensure_double_maths_files_clean() {
	local rel_file repo_file

	for rel_file in "${DOUBLE_MATHS_FILES[@]}"; do
		repo_file="src/KOKKOS/${rel_file#./}"
		if ! git -C "$LMPDIR" diff --quiet -- "$repo_file"; then
			echo "# Refusing double_maths_ops because $repo_file has local changes" | tee -a "$EXDIR/test_KOKKOS.log"
			return 1
		fi
	done
}

enable_double_maths_ops() {
	local rel_file

	if [[ "$USE_DOUBLE_MATHS_OPS" -ne 1 ]]; then
		return
	fi

	echo '#' | tee -a "$EXDIR/test_KOKKOS.log"
	echo '# double_maths_ops enabled: temporarily replacing sqrtf/expf with sqrt/exp' | tee -a "$EXDIR/test_KOKKOS.log"

	find_double_maths_files

	if [[ ${#DOUBLE_MATHS_FILES[@]} -eq 0 ]]; then
		echo '# No src/KOKKOS files currently contain sqrtf/expf. Nothing to override.' | tee -a "$EXDIR/test_KOKKOS.log"
		return
	fi

	ensure_double_maths_files_clean || exit 1

	trap 'restore_double_maths_state "$?"' EXIT
	RESTORE_DOUBLE_MATHS_OPS=1

	cd "$LMPDIR/src/KOKKOS" || exit 1
	echo '# The following files will be restored with git restore after the run:' | tee -a "$EXDIR/test_KOKKOS.log"
	for rel_file in "${DOUBLE_MATHS_FILES[@]}"; do
		echo "#   - $rel_file" | tee -a "$EXDIR/test_KOKKOS.log"
		sed -i -e 's/sqrtf(/sqrt(/g' -e 's/expf(/exp(/g' "$rel_file"
	done
	echo "# Total files modified: ${#DOUBLE_MATHS_FILES[@]}" | tee -a "$EXDIR/test_KOKKOS.log"
}

restore_double_maths_state() {
	local exit_code="$1"

	trap - EXIT
	restore_double_maths_now
	exit "$exit_code"
}

restore_double_maths_now() {
	local rel_file repo_file

	if [[ "$RESTORE_DOUBLE_MATHS_OPS" -ne 1 ]]; then
		return
	fi

	echo '# Restoring original sqrtf/expf calls with git restore' | tee -a "$EXDIR/test_KOKKOS.log"

	for rel_file in "${DOUBLE_MATHS_FILES[@]}"; do
		repo_file="src/KOKKOS/${rel_file#./}"
		git -C "$LMPDIR" restore -- "$repo_file" || return 1
	done

	RESTORE_DOUBLE_MATHS_OPS=0
}

###################################################################################################
# NVE/oxdna3-NVT Helper Functions
###################################################################################################
compare_nve_data() {
	local ref_file="$1"
	local test_file="$2"
	local tol="$3"
	local label="$4"
	local mode="$5"

	paste "$ref_file" "$test_file" |
		awk -v tol="$tol" -v label="$label" -v mode="$mode" '
			function check_pair(i, j,    diff) {
				diff = 0
				if ($(i) != 0) diff = ($(i)-$(j))/$(i)
				if (diff < 0) diff = -diff
				if (diff > tol) {
					printf "# Line %d: %g vs %g (relative difference = %g > %g)\n", NR, $(i), $(j), diff, tol
					printf "# %s FAILED\n", label
					failed = 1
					exit 1
				}
			}
			failed == 0 {
				if (mode == "extended") {
					for (i = 4; i <= 48; i += 4) check_pair(i, i + 48)
				} else {
					check_pair(4, 20)
					check_pair(8, 24)
					check_pair(12, 28)
					check_pair(16, 32)
				}
			}
			END {
				if (failed == 0) printf "# %s passed\n", label
			}
		'
}

run_nve_backend_checks() {
	local input_file="$1"
	local log_stem="$2"
	local mode="$3"

	### 1 MPI-task ###
	mpirun -np 1 ./lmp -in "$input_file" -k on -sf kk -pk kokkos comm device > /dev/null
	mv log.lammps "log.$DATE.$log_stem.g++.1"
	grep -e '[0-9]  ekin' "log.$DATE.$log_stem.g++.1" > e_test.1.dat
	grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat
	compare_nve_data e_ref.1.dat e_test.1.dat "$REL_TOL" "1 MPI-task" "$mode" 2>&1 | tee -a "$EXDIR/test_KOKKOS.log"

	### 4 MPI-tasks ###
	mpirun -np 4 ./lmp -in "$input_file" -k on -sf kk -pk kokkos comm device > /dev/null
	mv log.lammps "log.$DATE.$log_stem.g++.4"
	grep -e '[0-9]  ekin' "log.$DATE.$log_stem.g++.4" > e_test.4.dat
	grep -e '[0-9]  ekin' ../log*4 > e_ref.4.dat
	compare_nve_data e_ref.4.dat e_test.4.dat "$REL_TOL" "4 MPI-tasks" "$mode" 2>&1 | tee -a "$EXDIR/test_KOKKOS.log"

	if [ "$HIP_TEST_FLAG" = true ]; then
		### HIP ###
		rm -rf ./lmp
		cp "$BUILDDIR_KK_HIP_OMP/lmp" .
		mpirun -np 1 ./lmp -in "$input_file" -k on g 1 t 2 -sf kk -pk kokkos neigh half > /dev/null
		mv log.lammps "log.$DATE.$log_stem.g++.hip.g1t2"
		grep -e '[0-9]  ekin' "log.$DATE.$log_stem.g++.hip.g1t2" > e_test.hip.g1t2.dat
		grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat
		compare_nve_data e_ref.1.dat e_test.hip.g1t2.dat "$REL_TOL_GPU" "HIP g1t2" "$mode" 2>&1 | tee -a "$EXDIR/test_KOKKOS.log"
	else
		echo '# Skipping HIP test block (HIP_TEST_FLAG=false)' | tee -a "$EXDIR/test_KOKKOS.log"
	fi

	### CUDA ###
	rm -rf ./lmp
	cp "$BUILDDIR_KK_CUDA_OMP/lmp" .
	mpirun -np 1 ./lmp -in "$input_file" -k on g 1 t 2 -sf kk -pk kokkos neigh half > /dev/null
	mv log.lammps "log.$DATE.$log_stem.g++.cuda.g1t2"
	grep -e '[0-9]  ekin' "log.$DATE.$log_stem.g++.cuda.g1t2" > e_test.cuda.g1t2.dat
	grep -e '[0-9]  ekin' ../log*1 > e_ref.1.dat
	compare_nve_data e_ref.1.dat e_test.cuda.g1t2.dat "$REL_TOL_GPU" "CUDA g1t2" "$mode" 2>&1 | tee -a "$EXDIR/test_KOKKOS.log"
}

run_nve_case() {
	local title="$1"
	local case_dir="$2"
	local input_file="$3"
	local data_file="$4"
	local log_stem="$5"
	local mode="$6"
	local fixed_extra_file="${7:-}"
	local lj_extra_file="${8:-}"
	local real_extra_file="${9:-}"

	printf '\n# Running %s\n' "$title" | tee -a "$EXDIR/test_KOKKOS.log"
	cd "$EXDIR/$case_dir" || exit 1
	mkdir test
	cd test || exit 1
	cp "$BUILDDIR_KK_SERIAL/lmp" .
	cp "../$input_file" .
	cp "../$data_file" .

	if [ -n "$fixed_extra_file" ]; then
		cp "../$fixed_extra_file" .
	fi

	if [ -n "$lj_extra_file" ]; then
		if [ "$UNITS" = lj ]; then
			cp "../$lj_extra_file" .
		elif [ "$UNITS" = real ] && [ -n "$real_extra_file" ]; then
			cp "../$real_extra_file" .
		fi
	fi

	run_nve_backend_checks "$input_file" "$log_stem" "$mode"
}

run_oxdna3_nvt_case() {
	printf '\n# Running oxDNA3 NVT and unique base pairing test\n' | tee -a $EXDIR/test_KOKKOS.log
	cd $EXDIR/oxDNA3/unique_bp
	mkdir test
	cd test
	cp $BUILDDIR_KK_SERIAL/lmp .
	cp ../in.dsring2 .
	cp ../data.dsring2 .
	cp ../oxdna3_lj.cgdna .
	mpirun -np 8 ./lmp -in in.dsring2 -k on -sf kk -pk kokkos comm device > /dev/null
	mv log.lammps log.$DATE.dsring2.g++.8
	grep -e '[0-9]  ekin' log.$DATE.dsring2.g++.8 | awk '{print $1, $12}' > edyn_test.8.dat
	grep -e '[0-9]  ekin' log.$DATE.dsring2.g++.8 | awk '{print $1, $28}' > ehbond_test.8.dat
	grep -e '[0-9]  ekin' ../log*dsring2.g++.8 | awk '{print $1, $12}' > edyn_ref.8.dat
	grep -e '[0-9]  ekin' ../log*dsring2.g++.8 | awk '{print $1, $28}' > ehbond_ref.8.dat
	avg_edyn_test=$(awk '{sum += $2; n++} END {if (n > 0) print sum / n}' edyn_test.8.dat)
	avg_edyn_ref=$(awk '{sum += $2; n++} END {if (n > 0) print sum / n}' edyn_ref.8.dat)
	avg_ehbond_test=$(awk 'NR > 2000 {sum += $2; n++} END {if (n > 0) print sum / n}' ehbond_test.8.dat)
	avg_ehbond_ref=$(awk 'NR > 2000 {sum += $2; n++} END {if (n > 0) print sum / n}' ehbond_ref.8.dat)
	tol=$REL_TOL_NVT
	ekin=44.4
	diff=$(echo "($avg_edyn_test - $ekin)/$ekin" | bc -l)
	diff=$(echo "if ($diff < 0) -1 * $diff else $diff" | bc -l)
	if (( $(echo "$diff > $REL_TOL_NVT" | bc -l) )); then
		printf "# Relative difference of kinetic energy %g > %g\n" "$diff" "$tol" | tee -a $EXDIR/test_KOKKOS.log
		echo '# 8 MPI-tasks NVT FAILED' | tee -a $EXDIR/test_KOKKOS.log
	else
		echo '# 8 MPI-tasks NVT passed' | tee -a $EXDIR/test_KOKKOS.log
	fi
	diff=$(echo "($avg_ehbond_test - $avg_ehbond_ref)/$avg_ehbond_ref" | bc -l)
	diff=$(echo "if ($diff < 0) -1 * $diff else $diff" | bc -l)
	if (( $(echo "$diff > $REL_TOL_NVT" | bc -l) )); then
		printf "# Relative difference of hydrogen bonding energy %g > %g\n" "$diff" "$tol" | tee -a $EXDIR/test_KOKKOS.log
		echo '# 8 MPI-tasks unique base pairing FAILED' | tee -a $EXDIR/test_KOKKOS.log
	else
		echo '# 8 MPI-tasks unique base pairing passed' | tee -a $EXDIR/test_KOKKOS.log
	fi
  # GPU tests probably not needed here, since the GPU RNGs are totally different to the CPU reference.
}