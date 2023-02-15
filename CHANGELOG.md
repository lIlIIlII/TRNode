# CHANGELOG

## [3.1] - 2023-2-14

### Feature

* 1. add version control .git folder

### Fixed

* 1. Fix free operation in monitor function to avoid causing segment fault in the processor fetch behavior.
* 2. Fix incorrect calculation of progress and estimated time when encounter an unadaptable network.

### Refactored

* 1. Improve slot utilization by reducing malloc operation on different eta_qty.