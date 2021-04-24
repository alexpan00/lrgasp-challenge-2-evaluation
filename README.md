# encode_quantification
## Requirements
You should have python 3.6+ and linux installed on your machine
## Install
```
git clone https://github.com/Tidesun/encode_quantification.git
cd encode_quantification
python -m venv base
source base/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
## Usage
### Participant website
```
source base/bin/activate
python encode_quantification/index.py
```
Open the link in your browser after the prompt shown in your terminal \
Annotation results and true expression files should be in tsv format and annotation should be in gtf format.
### Organizer reporter
```
source base/bin/activate
cd encode_quantification
python -m report_generator.main -a ANNOTATION -r RESULT -t TRUTH -o OUTPUT --num_method NUM_METHOD  --num_samples NUM_SAMPLES
```
```
Quantification evaluation reporter

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -a ANNOTATION, --annotation ANNOTATION
                        The path of annotation file [GTF]
  -r RESULT, --result RESULT
                        The path of quantification result file [TSV\ZIP]
  -t TRUTH, --truth TRUTH
                        The path of true expression file [TSV]
  -o OUTPUT, --output OUTPUT
                        The path of output directory
  --num_method NUM_METHOD
                        Whether multi method data given ['Single' or 'Multi']
  --num_samples NUM_SAMPLES
                        Whether multi sample data given ['Single' or 'Multi']
```
