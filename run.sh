rm -rf /pnfs/sbnd/scratch/users/rlazur/MCP2025B/ccpi0/pot_5e18/batch
python3 batch/medulla.py -t selection/toml/grid_test.toml -p /pnfs/sbnd/scratch/users/rlazur/MCP2025B/ccpi0/pot_5e18/batch -b 1 --create-project
python3 batch/medulla.py -p /pnfs/sbnd/scratch/users/rlazur/MCP2025B/ccpi0/pot_5e18/batch --test-job
