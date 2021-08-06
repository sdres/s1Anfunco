#!/usr/bin/python
import os


def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes


def infotodict(seqinfo):



    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """

    rest3b = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-restBA3b_run-00{item:01d}_cbv')
    rest1 = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-restBA1_run-00{item:01d}_cbv')
    stim = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-stim_run-00{item:01d}_cbv')
    epi_stim = create_key('sub-{subject}/{session}/func/sub-{subject}_{session}_task-stim_run-00{item:01d}_bold')
    slab_inv1 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_inv1')
    slab_inv1_phs = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_inv1_phs')
    slab_inv2 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_inv2')
    slab_inv2_phs=create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_inv2_phs')
    slab_div=create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_div')
    slab_uni=create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_uni')
    slab_t1=create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_highres-mp2rage_run-00{item:01d}_t1')
    wholebrain_inv1 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_inv1')
    wholebrain_inv1_phs = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_inv1_phs')
    wholebrain_inv2 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_inv2')
    wholebrain_inv2_phs = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_inv2_phs')
    wholebrain_div = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_div')
    wholebrain_uni =create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_uni')
    wholebrain_t1 = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_t1')
    wholebrain_uni_den = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_wholebrain-mp2rage_run-00{item:01d}_uni-den')

    info = {rest3b:[],rest1:[],stim:[],epi_stim:[],slab_inv1:[],slab_inv1_phs:[],slab_inv2:[],slab_inv2_phs:[],slab_div:[],slab_uni:[],slab_t1:[],wholebrain_inv1:[], wholebrain_inv1_phs:[],wholebrain_inv2:[],wholebrain_inv2_phs:[],wholebrain_div:[],wholebrain_uni:[],wholebrain_t1:[],wholebrain_uni_den:[]}

    for idx, s in enumerate(seqinfo):
        if ('VASO_thick_slices_BA3b' in s.protocol_name) and (s.series_files == 318):
            info[rest3b].append(s.series_id)
        if ('VASO_thick_slices_BA1_rest_E00_M' in s.protocol_name):
            info[rest1].append(s.series_id)
        if ('VASO_thick_slices_BA3b_task' in s.protocol_name) and (s.series_files >= 200):
            info[stim].append(s.series_id)
        if ('EPI_wihout_inversion' in s.protocol_name) and (s.series_files >= 200):
            info[epi_stim].append(s.series_id)
        if ('High_Res_Slab_Anat_INV1' in s.series_description):
            info[slab_inv1].append(s.series_id)
        if ('High_Res_Slab_Anat_INV1_PHS' in s.series_description):
            info[slab_inv1_phs].append(s.series_id)
        if ('High_Res_Slab_Anat_INV2' in s.series_description):
            info[slab_inv2].append(s.series_id)
        if ('High_Res_Slab_Anat_INV2_PHS' in s.series_description):
            info[slab_inv2_phs].append(s.series_id)
        if ('High_Res_Slab_Anat_DIV_Images' in s.series_description):
            info[slab_div].append(s.series_id)
        if ('High_Res_Slab_Anat_UNI_Images' in s.series_description):
            info[slab_uni].append(s.series_id)
        if ('MP2RAGE_INV1' in s.series_description):
            info[wholebrain_inv1].append(s.series_id)
        if ('MP2RAGE_INV1_PHS' in s.series_description):
            info[wholebrain_inv1_phs].append(s.series_id)
        if ('MP2RAGE_INV2' in s.series_description):
            info[wholebrain_inv2].append(s.series_id)
        if ('MP2RAGE_INV2_PHS' in s.series_description):
            info[wholebrain_inv2_phs].append(s.series_id)
        if ('MP2RAGE_DIV_Images' in s.series_description):
            info[wholebrain_div].append(s.series_id)
        if ('MP2RAGE_UNI_Images' in s.series_description):
            info[wholebrain_uni].append(s.series_id)
        if ('MP2RAGE_T1_Images' in s.series_description):
            info[wholebrain_t1].append(s.series_id)
        if ('MP2RAGE_UNI-DEN' in s.series_description):
            info[wholebrain_uni_den].append(s.series_id)
    return info
