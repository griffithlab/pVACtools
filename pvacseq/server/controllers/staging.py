import os
import subprocess
import json
from flask import current_app
from yaml import dump
from .utils import initialize, savedata


def staging(input, samplename, alleles, epitope_lengths, prediction_algorithms,
          peptide_sequence_length, gene_expn_file, transcript_expn_file,
          normal_snvs_coverage_file, normal_indels_coverage_file,
          tdna_snvs_coverage_file, tdna_indels_coverage_file,
          trna_snvs_coverage_file, trna_indels_coverage_file,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold, fasta_size, iedb_retries, downstream_sequence_length, keep_tmp_files):
    """Stage input for a new pVAC-Seq run.  Generate a unique output directory and \
    save uploaded files to temporary locations (and give pVAC-Seq the filepaths). \
    Then forward the command to start()"""
    data = initialize()
    current_path = os.path.join(os.path.expanduser('~'), "Documents", "pVAC-Seq Output", samplename)
    if os.path.exists(current_path):
        i = 1
        while os.path.exists(current_path+"_"+str(i)):
            i += 1
        current_path += "_"+str(i)

    os.makedirs(os.path.join(current_path, 'Staging'), exist_ok=True)

    staged_input = open(os.path.join(current_path, "Staging", "input.vcf"), 'wb')
    input.save(staged_input)
    staged_input.flush()

    staged_additional_input_file_list = open(os.path.join(current_path, "Staging", "additional_input_file_list.yml"), 'w')

    staged_gene_expn_file = open(os.path.join(current_path, "Staging", "genes.fpkm_tracking"), 'wb')
    gene_expn_file.save(staged_gene_expn_file)
    staged_gene_expn_file.flush()
    if staged_gene_expn_file.tell():
        dump({"gene_expn_file": staged_gene_expn_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_transcript_expn_file = open(os.path.join(current_path, "Staging", "transcript.fpkm_tracking"), 'wb')
    transcript_expn_file.save(staged_transcript_expn_file)
    staged_transcript_expn_file.flush()
    if staged_transcript_expn_file.tell():
        dump({"transcript_expn_file" : staged_transcript_expn_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_normal_snvs_coverage_file = open(os.path.join(current_path, "Staging", "normal_snvs.bam_readcount"), 'wb')
    normal_snvs_coverage_file.save(staged_normal_snvs_coverage_file)
    staged_normal_snvs_coverage_file.flush()
    if staged_normal_snvs_coverage_file.tell():
        dump({"normal_snvs_coverage_file" : staged_normal_snvs_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_normal_indels_coverage_file = open(os.path.join(current_path, "Staging", "normal_indels.bam_readcount"), 'wb')
    normal_indels_coverage_file.save(staged_normal_indels_coverage_file)
    staged_normal_indels_coverage_file.flush()
    if staged_normal_indels_coverage_file.tell():
        dump({"normal_indels_coverage_file" : staged_normal_indels_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_tdna_snvs_coverage_file = open(os.path.join(current_path, "Staging", "tdna_snvs.bam_readcount"), 'wb')
    tdna_snvs_coverage_file.save(staged_tdna_snvs_coverage_file)
    staged_tdna_snvs_coverage_file.flush()
    if staged_tdna_snvs_coverage_file.tell():
        dump({"tdna_snvs_coverage_file" : staged_tdna_snvs_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_tdna_indels_coverage_file = open(os.path.join(current_path, "Staging", "tdna_indels.bam_readcount"), 'wb')
    tdna_indels_coverage_file.save(staged_tdna_indels_coverage_file)
    staged_tdna_indels_coverage_file.flush()
    if staged_tdna_indels_coverage_file.tell():
        dump({"tdna_indels_coverage_file" : staged_tdna_indels_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_trna_snvs_coverage_file = open(os.path.join(current_path, "Staging", "trna_snvs.bam_readcount"), 'wb')
    trna_snvs_coverage_file.save(staged_trna_snvs_coverage_file)
    staged_trna_snvs_coverage_file.flush()
    if staged_trna_snvs_coverage_file.tell():
        dump({"trna_snvs_coverage_file" : staged_trna_snvs_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_trna_indels_coverage_file = open(os.path.join(current_path, "Staging", "trna_indels.bam_readcount"), 'wb')
    trna_indels_coverage_file.save(staged_trna_indels_coverage_file)
    staged_trna_indels_coverage_file.flush()
    if staged_trna_indels_coverage_file.tell():
        dump({"trna_indels_coverage_file" : staged_trna_indels_coverage_file.name}, staged_additional_input_file_list)
    staged_additional_input_file_list.flush()

    return start(staged_input.name, samplename, alleles, epitope_lengths, prediction_algorithms, current_path,
              peptide_sequence_length, staged_additional_input_file_list.name if staged_additional_input_file_list.tell() else "", # check if any data written to file
              net_chop_method, len(netmhc_stab), len(top_result_per_mutation), top_score_metric,
              binding_threshold, minimum_fold_change,
              normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
              expn_val, net_chop_threshold,
              fasta_size, iedb_retries, downstream_sequence_length, len(keep_tmp_files))


def start(input, samplename, alleles, epitope_lengths, prediction_algorithms, output,
          peptide_sequence_length, additional_input_file_list,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold,
          fasta_size, iedb_retries, downstream_sequence_length, keep_tmp_files):
    """Build the command for pVAC-Seq, then spawn a new process to run it"""
    data = initialize()

    command = [
        'pvacseq',
        'run',
        input,
        samplename,
        alleles
    ]
    for algo in prediction_algorithms.split(','):
        command.append(algo)
    command += [
        output,
        '-e', epitope_lengths,
        '-l', str(peptide_sequence_length),
        '-m', top_score_metric,
        '-b', str(binding_threshold),
        '-c', str(minimum_fold_change),
        '--normal-cov', str(normal_cov),
        '--tdna-cov', str(tdna_cov),
        '--trna-cov', str(trna_cov),
        '--normal-vaf', str(normal_vaf),
        '--tdna-vaf', str(tdna_vaf),
        '--trna-vaf', str(trna_vaf),
        '--expn-val', str(expn_val),
        '-s', str(fasta_size),
        '-r', str(iedb_retries),
        '-d', str(downstream_sequence_length)
    ]
    if len(additional_input_file_list):
        command += ['-i', additional_input_file_list]
    if len(net_chop_method):
        command += [
            '--net-chop-method', net_chop_method,
            '--net-chop-threshold', str(net_chop_threshold)
        ]
    if netmhc_stab:
        command.append('--netmhc-stab')
    if top_result_per_mutation:
        command.append('--top-result-per-mutation')
    if keep_tmp_files:
        command.append('-k')

    # stdout and stderr from the child process will be directed to this file
    logfile = os.path.join(output, 'pVAC-Seq.log')

    data['processid']+=1
    os.makedirs(os.path.dirname(logfile), exist_ok = True)
    current_app.config['storage']['children'][data['processid']] = subprocess.Popen(
        command,
        stdout=open(logfile, 'w'),  # capture stdout in the logfile
        stderr=subprocess.STDOUT,
        # isolate the child in a new process group
        # this way it will remainin running no matter what happens to this process
        preexec_fn=os.setpgrp
    )
    # Store some data about the child process
    data['process-%d'%(data['processid'])] = {
        # Do the replacement so that the displayed command is actually valid
        # The executed command is automatically escaped as part of Popen
        'command': " ".join(command),
        'logfile':logfile,
        'pid':current_app.config['storage']['children'][data['processid']].pid,
        'status': "Task Started",
        'output':os.path.abspath(output)
    }
    if 'reboot' not in data:
        data['reboot'] = current_app.config['reboot']
    savedata(data)
    configObj = {
        'action':'run',
        'input_file': input,
        'sample_name':samplename,
        'alleles':','.split(alleles),
        'prediction_algorithms':','.split(prediction_algorithms),
        'output_directory':output
    }
    if epitope_lengths!=10:
        configObj['epitope_lengths']=','.split(epitope_lengths)
    if peptide_sequence_length!=21:
        configObj['peptide_sequence_length']=peptide_sequence_length
    if additional_input_file_list!='':
        configObj['additional_input_files']=','.split(additional_input_file_list)
    if net_chop_method!='':
        configObj['net_chop_method']=net_chop_method
    if netmhc_stab:
        configObj['netmhc_stab']=True
    if top_result_per_mutation:
        configObj['top_result_per_mutation']=True
    if top_score_metric!='median':
        configObj['top_score_metric']=top_score_metric
    if binding_threshold!=500:
        configObj['binding_threshold']=binding_threshold
    if minimum_fold_change!=0:
        configObj['minimum_fold_change']=minimum_fold_change
    if normal_cov!=5:
        configObj['normal_coverage_cutoff']=normal_cov
    if tdna_cov!=10:
        configObj['tumor_dna_coverage_cutoff']=tdna_cov
    if trna_cov!=10:
        configObj['tumor_rna_coverage_cutoff']=trna_cov
    if normal_vaf!=2:
        configObj['normal_vaf_cutoff']=normal_vaf
    if tdna_vaf!=40:
        configObj['tumor_dna_vaf_cutoff']=tdna_vaf
    if trna_vaf!=40:
        configObj['tumor_rna_vaf_cutoff']=trna_vaf
    if expn_val!=1:
        configObj['expression_cutoff']=expn_val
    if net_chop_threshold!=0.5:
        configObj['netchop_threshold']=net_chop_threshold
    if fasta_size!=200:
        configObj['fasta_size']=fasta_size
    if iedb_retries!=5:
        configObj['iedb_retries'] = iedb_retries
    if downstream_sequence_length!="1000":
        configObj['downstream_sequence_length']=downstream_sequence_length

    writer = open(os.path.join(
        os.path.abspath(output),
        'config.json'
    ),'w')
    json.dump(configObj, writer, indent='\t')
    writer.close()
    return data['processid']


def test():
    """Return the submission page (a stand-in until there is a proper ui for submission)"""
    reader = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_start.html'))
    data = reader.read()
    reader.close()
    return data


def check_allele(allele):
    """Checks if the requested allele is supported by pVAC-Seq or not"""
    data = initialize()
    if 'allele_file' not in current_app.config:
        allele_file = tempfile.TemporaryFile('w+')
        subprocess.call(['pvacseq', 'valid_alleles'], stdout=allele_file)
        current_app.config['allele_file'] = allele_file
    current_app.config['allele_file'].seek(0)
    for line in current_app.config['allele_file']:
        if line.strip() == allele:
            return True
    return False
