for ligand in ligands/*.pdbqt
do
    name=$(basename "$ligand" .pdbqt)

    echo "======================================="
    echo "Docking: $name"
    echo "======================================="

    vina \
      --config config.txt \
      --ligand "$ligand" \
      --out "docking_output/${name}_out.pdbqt" \
      > "docking_logs/${name}.log"

    echo "Finished: $name"
done

echo "======================================="
echo "ALL DOCKING COMPLETED!"
echo "======================================="
