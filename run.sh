rm -rf /pnfs/sbnd/scratch/users/rlazur/MCP2025B/ccpi0/pot_5e18/batch
python3 batch/medulla.py --toml selection/toml/grid_test.toml --project-dir /pnfs/sbnd/scratch/users/rlazur/MCP2025B/ccpi0/pot_5e18/batch --batch-size 1 --create-project
python3 batch/medulla.py --project-dir /pnfs/sbnd/scratch/users/rlazur/MCP2025B/ccpi0/pot_5e18/batch --experiment sbnd --test-job
