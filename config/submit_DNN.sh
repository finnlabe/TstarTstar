#!/bin/bash

sframe_batch.py -s DNN_TstarTstar_2016.xml
sframe_batch.py -s DNN_TstarTstar_2017.xml
sframe_batch.py -s DNN_TstarTstar_2018.xml

sframe_batch.py -s DNN_SM_2016.xml
sframe_batch.py -s DNN_SM_2017.xml
sframe_batch.py -s DNN_SM_2018.xml

sframe_batch.py -s DNN_addSM_2016.xml
sframe_batch.py -s DNN_addSM_2017.xml
sframe_batch.py -s DNN_addSM_2018.xml

sframe_batch.py -s DNN_Data_2016.xml
sframe_batch.py -s DNN_Data_2017.xml
sframe_batch.py -s DNN_Data_2018.xml

sframe_batch.py -s DNN_DataDrivenBG_2016.xml
sframe_batch.py -s DNN_DataDrivenBG_2017.xml
sframe_batch.py -s DNN_DataDrivenBG_2018.xml
