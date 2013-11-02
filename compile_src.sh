#!/bin/bash

R CMD SHLIB CCM_bootstrap.c
R CMD SHLIB SSR_predict_boot.c
R CMD SHLIB bmod.c
