#!/bin/bash

export OMP_NUM_THREADS=20

cd ~/ghastly/openmc_models/inf_lat_dep

cd iter0
python3 bcc_model.py
cp depletion_results.h5 ../passwise/iter1/i0-dep-res.h5
cp depletion_results.h5 ../corewise/iter1/i0-dep-res.h5

cd ../passwise
python3 gen_inputs.py

cd iter1
python3 bcc_model.py
cp depletion_results.h5 ../iter2/i1-dep-res.h5

cd ../iter2
python3 bcc_model.py
cp depletion_results.h5 ../iter3/i2-dep-res.h5

cd ../iter3
python3 bcc_model.py
cp depletion_results.h5 ../iter4/i3-dep-res.h5

cd ../iter4
python3 bcc_model.py
cp depletion_results.h5 ~/ghastly/openmc_models/bcc_phys/passwise/i4-dep-res.h5

cd ~/ghastly/openmc_models/inf_lat_dep

cd ../corewise
python3 gen_inputs.py

cd iter1/
python3 bcc_model.py
cp depletion_results.h5 ../iter2/i1-dep-res.h5

cd ../iter2
python3 bcc_model.py
cp depletion_results.h5 ../iter3/i2-dep-res.h5

cd ../iter3
python3 bcc_model.py
cp depletion_results.h5 ../iter4/i3-dep-res.h5

cd ../iter4
python3 bcc_model.py
cp depletion_results.h5 ~/ghastly/openmc_models/bcc_phys/corewise/i4-dep-res.h5

cd ~/ghastly/openmc_models/bcc_phys

cd passwise
python3 gen_inputs.py

cd /0days/
python3 bcc_model.py

cd ../249days/
python3 bcc_model.py

cd ../499days/
python3 bcc_model.py

cd ../749days/
python3 bcc_model.py

cd ../999days/
python3 bcc_model.py

cd ../1249days/
python3 bcc_model.py

cd ../1549days/
python3 bcc_model.py

cd ~/ghastly/openmc_models/bcc_phys

cd corewise
python3 gen_inputs.py

cd /0days/
python3 bcc_model.py

cd ../249days/
python3 bcc_model.py

cd ../499days/
python3 bcc_model.py

cd ../749days/
python3 bcc_model.py

cd ../999days/
python3 bcc_model.py

cd ../1249days/
python3 bcc_model.py

cd ../1549days/
python3 bcc_model.py


