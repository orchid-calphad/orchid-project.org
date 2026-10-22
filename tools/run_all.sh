#!/bin/zsh
# Full SGTE binary phase-diagram build (2026-10-22): diagrams -> grids -> sources -> index
set -u
PY=~/orchid-py313/bin/python; SITE=~/projects/orchid-project.org; T=~/projects/open-tdb/databases/sgte-binary
OUT=$SITE/data/phase-diagrams/sgte; LOG=$SITE/tools/run_all.log
cd $T
echo "start $(date)" > $LOG
$PY $SITE/tools/compute_binary_diagrams.py $OUT $(ls | grep -i "\.tdb$") --workers 9 2>/dev/null | grep -E " (ok|failed|exists) " >> $LOG
echo "diagrams done $(date)" >> $LOG
$PY $SITE/tools/add_phase_grid.py $OUT $T --workers 9 2>/dev/null | grep -E " (ok|failed|skip) " >> $LOG
echo "grids done $(date)" >> $LOG
python3 $SITE/tools/add_sources_and_index.py >> $LOG 2>&1
echo "all done $(date)" >> $LOG
