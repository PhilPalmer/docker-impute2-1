params.rawdata               = "${baseDir}/data/id_600_raw_data.zip"
params.chromosome            = (21..22) + 'X'
params.uniqueID              = "placeholder"
params.email                 = "alaincoletta@insilicodb.com"
params.protect_from_deletion = "TRUE"
params.filename              = "placeholder.txt"
params.webhook               = false
params.isFileValidation      = false
params.publishDir            = "/tmp"
params.referenceUrl          = "http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz"
params.resultdir             = "results"
uniqueID = params.uniqueID

if ( ! params.chromosome instanceof List) {
  params.chromosome = params.chromosome.tokenize(',')
}

if (uniqueID == "placeholder") {
  uniqueID = (params.rawdata =~ /([^\/]*)_raw_data.[a-z]+$/)[ 0 ][ 1 ]
}
println uniqueID
println params.rawdata

rawFileChannel = Channel.fromPath(params.rawdata)

process fetch_reference {
  echo true

  input:
  val referenceUrl from params.referenceUrl

  output:
  file('ALL_1000G_phase1integrated_v3_impute') into referenceDir

  script:
  """
  axel --quiet http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz

  tar xf ALL_1000G_phase1integrated_v3_impute.tgz
  rm ALL_1000G_phase1integrated_v3_impute.tgz
  cd ALL_1000G_phase1integrated_v3_impute

  mv genetic_map_chrX_nonPAR_combined_b37.txt genetic_map_chrX_combined_b37.txt
  mv ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz ALL_1000G_phase1integrated_v3_chrX_impute.hap.gz
  mv ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz ALL_1000G_phase1integrated_v3_chrX_impute.legend.gz
  """

}

referenceDir1 = Channel.create()
referenceDir2 = Channel.create()
(referenceDir1, referenceDir2) = referenceDir.into(2)

process unzip {
  echo true

  input:
  file inputFile from rawFileChannel

  output:
  file("*.txt") into unzipOutChan

  script:
  if ( "${inputFile}".endsWith(".zip") ) {
    """
    mkdir zip_extract
    unzip ${inputFile} -d zip_extract
    pushd zip_extract
    for file in `ls -p | grep -v / `;
    do
      if [ \${file: -4} == ".txt" ] || [ \${file: -4} == ".csv" ];
      then
        echo \$file
        echo cp \$file ../rawdata.txt
        cp \$file ../rawdata.txt
      fi
    done
    popd
    """
  } else if ( "${inputFile}".endsWith(".txt") ) {
    """
    cp ${inputFile} rawdata.txt
    """
  } else if ( "${inputFile}".endsWith(".gz") ) {
    """
    gunzip -c ${inputFile} > rawdata.txt
    """
  } else {
    error "Unsupported file extension ${inputFile}"
  }
}

process convert {
  echo true

  input:
  file inputFile from unzipOutChan

  output:
  file("sorted_rawData.txt") into convertChan
  file("sorted_rawData.txt") into convertChan2
  file("X_rawData.txt") into convertChan3

  """
  # Remove #, rsid
  # Replace comma with tabs
  # Replace quotes
  # Replace chromosome 23 with X, 24 with Y, 25 with XY, 26 with MT
  # Concat column 4 and 5 in case of Ancestry.dna
  grep -v "#" ${inputFile} |
  sed "s/,/\\t/g" |
  sed "s/\\"//g" |
  grep -Piv "rsid\\t" |
  sed -E "s/([a-z]+[0-9]+)\\t23\\t/\\1\\tX\\t/g" |
  sed -E "s/([a-z]+[0-9]+)\\t24\\t/\\1\\tY\\t/g" |
  sed -E "s/([a-z]+[0-9]+)\\t25\\t/\\1\\tXY\\t/g" |
  sed -E "s/([a-z]+[0-9]+)\\t26\\t/\\1\\tMT\\t/g" |
  awk '{ print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 \$5}' > tempFile.txt

  # sort content of file on chromosome, needed for plink
  # disable errors for headers and X,Y, MT because doesn't need to be there
  set +e
  grep "^#" ${inputFile} > sorted_rawData.txt
  set -e
  for i in {1..22}; do grep -P "^[a-z]+[0-9]+\\t\$i\\t" tempFile.txt >> sorted_rawData.txt; done

  set +e
  grep -P "^[a-z]+[0-9]+\\tX\\t" tempFile.txt > X_rawData.txt
  if [[ -s X_rawData.txt ]]; then
    cat X_rawData.txt >> sorted_rawData.txt;
  else
    cat /scripts/data/X_dummy_data.txt >> sorted_rawData.txt;
  fi
  grep -P "^[a-z]+[0-9]+\\tY\\t" tempFile.txt >> sorted_rawData.txt
  grep -P "^[a-z]+[0-9]+\\tXY\\t" tempFile.txt >> sorted_rawData.txt
  grep -P "^[a-z]+[0-9]+\\tMT\\t" tempFile.txt >> sorted_rawData.txt
  set -e
  """
}

process plink {

  echo true
  maxForks 6

  input:
  file rawdata from convertChan

  output:
  file("step_1.{map,ped,log,hh}") into plinkOutChan
  file("step_1.{map,ped,log,hh}") into plinkOutChan2
  file("step_1.{map,ped,log,hh}") into plinkOutChan3
  file('step_2_exclusions') into exclusionsChan
  file('step_2_exclusions') into exclusionsChan2

  """
  plink-1.9 --noweb --23file  ${rawdata} John Doe --recode --out step_1
  Rscript '/scripts/remove-duplicates.R' step_1.map step_2_exclusions
  """
}

process createFileValidationOutput {

  echo true
  publishDir params.publishDir, mode: 'copy'

  when:
  params.isFileValidation

  input:
  file("*") from plinkOutChan3
  file('step_2_exclusions') from exclusionsChan2
  file('X_rawData.txt') from convertChan3

  output:
  file("validated.txt") into createFileValidationOutputChan

  """
  if [[ -s  X_rawData.txt ]]; then
    echo "valid" > validated.txt
  else
    echo "missing-x-chr" > validated.txt
  fi

  """
}

//duplicatesOutChan.subscribe{  println "${it}"}

process extractChromosome {

  echo true

  when:
  !params.isFileValidation

  input:
  each chromosome from ( params.chromosome )
  file("*") from plinkOutChan
  file('step_2_exclusions') from exclusionsChan

  output:
  set val(chromosome), file("step_2_chr*") into duplicatesOutChan

  """
  plink-1.9 --file step_1 --chr ${chromosome} --recode --out "step_2_chr${chromosome}" --exclude step_2_exclusions
  """
}

//extractChromosomeChan.subscribe{  println "${it}"}


lala1 = Channel.create()
(lala1, duplicatesRefOutChan) = referenceDir1.combine(duplicatesOutChan).into(2)

lala1.println()
//referenceDir1.combine(duplicatesOutChan).set { duplicatesRefOutChan } 

//duplicatesOutChan.map{ chromosome, step_files -> [ chromosome, file(params.referenceDir), step_files ] }.set { duplicatesRefOutChan } 

process checkStrandFlips {

  validExitStatus 0,1,2
  errorStrategy 'ignore'
  echo true

  input:
  set file('reference'), val(chr), file("*") from duplicatesRefOutChan

  output:
  set val(chr), file('maxLenght_chr*.txt'), file('step_5_chr*.haps'), file("step_4_chr${chr}.sample") into chromosomeMaxLengthOut

  script:
  genomeFile="reference/genetic_map_chr${chr}_combined_b37.txt"
  hapFile="reference/ALL_1000G_phase1integrated_v3_chr${chr}_impute.hap.gz"
  legendFile="reference/ALL_1000G_phase1integrated_v3_chr${chr}_impute.legend.gz"
  sampleFile="reference/ALL_1000G_phase1integrated_v3.sample"

  """
  echo "################### DOING ${chr} ####################"
  echo "################### DOING ${chr} ls BEGIN ####################"
  ls -l
  echo "################### DOING ${chr} ls END ####################"
  echo "shapeit -check --input-ped step_2_chr${chr}.ped step_2_chr${chr}.map -M $genomeFile --input-ref $hapFile $legendFile $sampleFile --output-log step_2_chr${chr}_shapeit_log"
  set +e
  shapeit -check --input-ped step_2_chr${chr}.ped step_2_chr${chr}.map -M $genomeFile --input-ref $hapFile $legendFile $sampleFile --output-log step_2_chr${chr}_shapeit_log
  set -e

  Rscript '/scripts/heterozygote.R' step_2_chr${chr}_shapeit_log.snp.strand step_2_chr${chr}.map step_2_chr${chr}.ped step_3_chr${chr}.ped step_3_chr${chr}_exclusions

  echo "shapeit --input-ped step_3_chr${chr}.ped step_2_chr${chr}.map -M $genomeFile --input-ref $hapFile $legendFile $sampleFile --output-log step_4_chr${chr}_shapeit_log --exclude-snp step_3_chr${chr}_exclusions -O step_4_chr${chr}"
  shapeit --input-ped step_3_chr${chr}.ped step_2_chr${chr}.map -M $genomeFile --input-ref $hapFile $legendFile $sampleFile --output-log step_4_chr${chr}_shapeit_log --exclude-snp step_3_chr${chr}_exclusions -O step_4_chr${chr}

  cut --delimiter=' ' -f 1-7 step_4_chr${chr}.haps > step_5_chr${chr}.haps

  head -n 3 step_4_chr${chr}.sample > step_5_chr${chr}.sample

  zcat $legendFile | tail -n 1 | cut --delimiter=\\  -f 2 > maxLenght_chr${chr}.txt

  ls -l
  """
}

chromosomeMaxLengthOut.flatMap(this.&foo).set { flattendChan }

referenceDir2.combine(flattendChan).set { flattendRefChan } 

//flattendChan.map{ chr, start, i, haps, sampleFile -> [ chr, file(params.referenceDir), start, i, haps, sampleFile ] }.set { flattendRefChan }

process impute2 {

  echo true
  validExitStatus 0,1,2
  errorStrategy 'ignore'
  maxForks 4

  input:
  set file('reference'), val(chr), val(start), val(i), file(haps), file(sampleFile) from flattendRefChan

  output:
  set val(chr), file("step_7_chr*"), file(sampleFile) into impute2Chan

  script:
  end = start + 5_000_000
  genomeFile="reference/genetic_map_chr${chr}_combined_b37.txt"
  hapFile="reference/ALL_1000G_phase1integrated_v3_chr${chr}_impute.hap.gz"
  legendFile="reference/ALL_1000G_phase1integrated_v3_chr${chr}_impute.legend.gz"

  """
  echo "impute2 -m $genomeFile -h $hapFile -l $legendFile -known_haps_g $haps -int ${start} ${end} -Ne 20000 -o step_7_chr${chr}_${i}"
  impute2 -m $genomeFile -h $hapFile -l $legendFile -known_haps_g $haps -int ${start} ${end} -Ne 20000 -o step_7_chr${chr}_${i}
  """

}

impute2GroupedChan = impute2Chan
  .flatMap{chr, files, sampleFile -> files.collect{[chr,it,sampleFile]}}
  .groupTuple()
  .groupTuple()
  .map {chr, files, sampleFile -> [chr, files.flatten().sort(), sampleFile]}
  .map {chr, files, sampleFile -> [chr, files.findAll { it=~/[0-9]$/ }, sampleFile]}
  .map {chr, files, sampleFile -> [chr, files.sort{(it.name=~ /.*_(\d+)$/)[ 0 ][ 1 ].toInteger()}, sampleFile.flatten().first()] }

//  .subscribe{  println "${it}"}


  /*
  .map {chr, files -> [chr, files.findAll { it=~/[0-9]$/ }, files.find {it=~/sample$/}]}
  .filter{chr, file -> file=~/[0-9]$/}
  .collectFile(sort:true)
  .map{file -> [file.name, file.renameTo("${uniqueID}_chr${file.name}.gen")]}
  */

process summarize {

  echo true

  input:
  set val(chr), file(genFiles), file(sampleFile) from impute2GroupedChan

  output:
  set file("${uniqueID}_chr${chr}.23andme.txt"), file("${uniqueID}_chr${chr}.gen") into summarizeChan

  script:
  genFile="${uniqueID}_chr${chr}.gen"

  """
  echo $chr
  echo $sampleFile
  echo $genFiles

  cat ${genFiles} > $genFile

  #make list of non-indels
  awk -F' ' '{ if ((length(\$4) > 1 ) || (length(\$5) > 1 )) print \$2 }' ${uniqueID}_chr${chr}.gen  > step_8_chr${chr}_snps_to_exclude

  /install/gtool -S --g ${uniqueID}_chr${chr}.gen --s $sampleFile --exclusion step_8_chr${chr}_snps_to_exclude --og step_8_chr${chr}.gen

  set +e
  #Convert to ped format
  /install/gtool -G --g step_8_chr${chr}.gen --s $sampleFile --chr ${chr} --snp

  #reform to plink fam/bim/bed file
  plink-1.07 --file step_8_chr${chr}.gen --recode --transpose --noweb --out step_9_chr${chr}

  awk '{ print \$2 \"\t\" \$1 \"\t\"\$4\"\t\" \$5 \$6}' step_9_chr${chr}.tped  > step_10_chr${chr}.txt
  set -e

  #The step 8 and also 9 sometime fails for no apparent reason. Probably memory. We therefore make a checkup, where
  #it is checked if the file actually exists and if not - a more complicated step splits it up in chunks.
  #It's not looking nice, but at least the split-up only needs to run in very low memory settings

  if [ -f step_10_chr${chr}.txt ]; then
    FILESIZE=\$(stat -c%s "step_10_chr${chr}.txt")
  else
    FILESIZE=0
  fi

  # arbitraly re-run if it's less than 100 bytes (fair to assume something was wrong then)
  if [ \$FILESIZE -lt 100]; then
   Rscript '/scripts/23andmeRerun.R' ${chr} "/install/gtool" "plink-1.07"
  fi

  #remove NN
  awk '{ if(\$4 != \"NN\") print}' step_10_chr${chr}.txt  > ${uniqueID}_chr${chr}.23andme.txt

  """

}
//summarizeChan.subscribe{it -> println it}

summarizeConcatChan = summarizeChan.reduce { a, b -> a + b }

process zipping {

  echo true

  when:
  !params.isFileValidation

  input:
  file toZipFiles from summarizeConcatChan
  file rawdata from convertChan2

  output:
  set file("${uniqueID}.gen.zip"), file("${uniqueID}.23andme.zip"), file("${uniqueID}.input_data.zip") into zippingOutChan

  script:
  zipFile23andme = "${uniqueID}.23andme.zip"
  zipFileGen = "${uniqueID}.gen.zip"
  // rawdata = rawdata.renameTo("${uniqueID}_raw_data.txt")

  """
  #zipping and moving 23andme files
  zip -r9X ${zipFile23andme} *.23andme.txt

  #zipping gen files
  zip -r9X ${zipFileGen} *.gen

  #move the original file as well
  zip -r9X ${uniqueID}.input_data.zip $rawdata
  """
}

// zippingOutChan.subscribe{it -> println it}

process createPDataFile {

  publishDir params.resultdir, mode: 'copy'

  echo true

  input:
  set file("${uniqueID}.gen.zip"), file("${uniqueID}.23andme.zip"), file("${uniqueID}.input_data.zip") from zippingOutChan
  file plinkFiles from plinkOutChan2

  output:
  set file("${uniqueID}.gen.zip"), file("${uniqueID}.23andme.zip"), file("${uniqueID}.input_data.zip"), file("pData.txt") into createPDataFileOutChan

  """
  Rscript '/scripts/createPDataFile.R' ${uniqueID} $params.email $params.filename $params.protect_from_deletion
  """
}

createPDataFileOutChan.subscribe{it -> println it}

workflow.onComplete {
  if (params.webhook && !workflow.success) {
    proc = ["curl", "-X", "POST", "-H", "Content-type: application/json", "--data", "'{\"text\": \"Imputation dataset ${uniqueID} failed\"}'", params.webhook].execute()
    proc.waitForProcessOutput(System.out, System.err)
  }
}


// --- helper methods

def foo( chr, maxFile, hapsFile, sampleFile ) {
  def time = System.currentTimeMillis()
  def max = maxFile.text.toBigInteger()

  def comb = [ [chr], range(0, max, 5_000_000) ].combinations()
  int index=1
  for ( row in comb ) {
    row << index++
    row << hapsFile
    row << sampleFile
  }
  return comb
}

def range(min, max, step) {
  def result = []
  def count = min;
  while( count<max ) {
      result.add(count)
      count+=step
  }
  return result
}


