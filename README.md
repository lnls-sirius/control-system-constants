# control-system-constants
Sirius Control System static parameters

This repository should be cloned as the 'html' directory at the Sirius control-system web server.

* beaglebone: static tables related to beaglebones.
* cs-studio-common: common CS-Studio files.
* diagnostics: data on diagnostics devices.
* documentation: general documentation on Sirius, including drawings, pictures, etc.
* macschedule: machine schedule files.
* magnet: static tables related to the magnets.
* pwrsupply: static tables related to magnet power supplies.
* timesys: static tables related to the timing subsystem.

## Deploy Procedure
The deploy of this webpage should be done in the IBM server 1 (IP address: 10.30.1.50) using docker compose.

IBM Server 1 Repository: https://gitlab.cnpem.br/SOL/Docker/server-ibm1/-/tree/master/fac-csconsts
