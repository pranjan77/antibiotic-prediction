# Antibiotic prediction

1. `train_classifiers.sh` - Script that trains the classifiers for BGC activity prediction (antibacterial, antifungal, etc.)
2. `predict_function.sh` - Script that predict the activity of BGCs for a given genome
3. `multiple_predict_function.py` - Script that runs antismash and rgi and predicts the activity of a list of genomes

Example command to run the workflow
```bash
python multiple_predict_function.py path/genome1 path/genome2 --output_dir outputs --no_SSN True
```
