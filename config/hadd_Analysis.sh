#!/bin/bash

sframe_batch.py -a Analysis_TstarTstar_2016.xml
sframe_batch.py -a Analysis_TstarTstar_2017.xml
sframe_batch.py -a Analysis_TstarTstar_2018.xml

sframe_batch.py -a Analysis_SM_2016.xml
sframe_batch.py -a Analysis_SM_2017.xml
sframe_batch.py -a Analysis_SM_2018.xml

sframe_batch.py -a Analysis_addSM_2016.xml
sframe_batch.py -a Analysis_addSM_2017.xml
sframe_batch.py -a Analysis_addSM_2018.xml

sframe_batch.py -a Analysis_Data_2016.xml
sframe_batch.py -a Analysis_Data_2017.xml
sframe_batch.py -a Analysis_Data_2018.xml
