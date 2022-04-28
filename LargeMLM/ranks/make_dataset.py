import os, random
import torch
import pandas as pd
from functools import partial
from datasets import load_dataset
from transformers import Trainer, TrainingArguments
from transformers import BigBirdModel
from transformers import BertTokenizerFast, PreTrainedTokenizerFast

from sklearn.model_selection import train_test_split

# Change this for other users:
#owens_desktop = '/Users/owenqueen/Desktop/bioinformatics/codonbert/CodonBert/Data'
#prefix = owens_desktop
prefix = '/lustre/isaac/scratch/oqueen/CodonBert/Data'

def prep_dataset():
    large_data = os.path.join(prefix, 'fna_txt')
    all_files = os.listdir(os.path.join(prefix, 'fna_txt'))
    all_files = [os.path.join(large_data, f) for f in all_files] # Transform to abs path

    val_pct = 0.1
    test_pct = 0.2
    batch_size = 256 # from RoBERTa
    seq_len = 512 # from RoBERTa

    # Choose val_pct of all_files for validation data:
    # val_sample = random.sample(range(len(all_files)), k = int(len(all_files) * val_pct))
    # val_files = [all_files[i] for i in val_sample]
    # val_set = set(val_sample) # Convert to set for O(1) lookup

    # # Choose test_pct of all_files for test data:
    # test_sample = random.sample(range(len(all_files)), k = int(len(all_files) * test_pct))
    # test_files = [all_files[i] for i in test_sample]
    # test_set = set(test_sample) # Convert to set for O(1) lookup

    train_inds, test_inds = train_test_split(list(range(len(all_files))), test_size = test_pct, random_state = None)
    train_inds, val_inds = train_test_split(train_inds, test_size = val_pct / (1 - test_pct), random_state = None)

    train_files = [all_files[i] for i in train_inds]
    test_files = [all_files[i] for i in test_inds]
    val_files = [all_files[i] for i in val_inds]

    file_dict = {
        'train': train_files,
        'validation': val_files,
        'test': test_files
    }

    dataset = load_dataset(
        'text',
        data_files = file_dict
    )

    print(dataset)
    print(type(dataset))

    # Load in wordlevel tokenizer:
    tokenizer = PreTrainedTokenizerFast(tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_rank.json'))

    def tokenize_function(examples):
        return tokenizer(examples['text'])

    tokenized_dataset = dataset.map(tokenize_function, batched = True, num_proc = 4, remove_columns=['text'])

    print(tokenized_dataset)
    print(tokenized_dataset['train'][0])

    def group_texts(examples, block_size = 128):
        # Concatenate all texts.
        concatenated_examples = {k: sum(examples[k], []) for k in examples.keys()}
        total_length = len(concatenated_examples[list(examples.keys())[0]])
        # We drop the small remainder, we could add padding if the model supported it instead of this drop, you can
            # customize this part to your needs.
        total_length = (total_length // block_size) * block_size
        # Split by chunks of max_len.
        result = {
            k: [t[i : i + block_size] for i in range(0, total_length, block_size)]
            for k, t in concatenated_examples.items()
        }
        result["labels"] = result["input_ids"].copy()
        return result

    gtexts = partial(group_texts, block_size = seq_len)

    lm_datasets = tokenized_dataset.map(gtexts, batched=True, batch_size = batch_size, num_proc=4)

    lm_datasets.save_to_disk('SavedDatasets/rank.hgf')

if __name__ == '__main__':
    prep_dataset()
