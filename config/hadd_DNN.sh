#!/bin/bash

sframe_batch.py -a DNN_TstarTstar_2016.xml
sframe_batch.py -a DNN_TstarTstar_2017.xml
sframe_batch.py -a DNN_TstarTstar_2018.xml

sframe_batch.py -a DNN_SM_2016.xml
sframe_batch.py -a DNN_SM_2017.xml
sframe_batch.py -a DNN_SM_2018.xml

sframe_batch.py -a DNN_addSM_2016.xml
sframe_batch.py -a DNN_addSM_2017.xml
sframe_batch.py -a DNN_addSM_2018.xml

sframe_batch.py -a DNN_Data_2016.xml
sframe_batch.py -a DNN_Data_2017.xml
sframe_batch.py -a DNN_Data_2018.xml
