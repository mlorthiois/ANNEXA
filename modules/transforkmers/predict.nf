process PREDICT {
  publishDir "$params.outdir/transforkmers", mode: 'copy'
  conda "pip pytorch::pytorch"
  cpus params.maxCpu
  memory '16GB'

  input:
  file tss_sequences

  output:
  path "output.csv", emit: tss_prob

  """
  pip install git+https://github.com/mlorthiois/transforkmers
  transforkmers predict \
    --model_path_or_name ${params.model} \
    --tokenizer ${params.tokenizer} \
    --input ${tss_sequences} \
    --quantize-model \
    --output . \
    --per_device_eval_batch_size 64
  """
}
