#BPH Trigger Tool 

A repository contains a set of "tests" that access the cms oms using cms api client and \
create a set of useful histograms. 


- Copy this repo : git clone 
- Install OMS API client : https://gitlab.cern.ch/cmsoms/oms-api-client. \
- Run the OMS API test :
  ```
  python3 oms-api-client/examples/01-query.py
  ```

- Run the BPH Trigger Tool:
  ```
  python3 my_test.py 
  ```
  This will run a default test : l1rate (l1 trigger rate per instant lumi/lumi section)

- For the full list of tests and flags run:
  ```
  python3 my_test.py --help
  ```


