import os, random
import torch
import pandas as pd
import numpy as np
from functools import partial
from datasets import load_dataset, load_from_disk
from transformers import Trainer, TrainingArguments
from transformers import BigBirdModel
from transformers import BertTokenizerFast, PreTrainedTokenizerFast

from sklearn.model_selection import train_test_split

# Change this for other users:
#owens_desktop = '/Users/owenqueen/Desktop/bioinformatics/codonbert/CodonBert/Data'
#prefix = owens_desktop
prefix = '/lustre/isaac/scratch/oqueen/CodonBert/Data'

val_pct = 0.1
test_pct = 0.2
batch_size = 256 # from RoBERTa
seq_len = 256 # from RoBERTa

def prep_dataset():
    large_data = os.path.join(prefix, 'fna_rankC')
    all_files = os.listdir(os.path.join(prefix, 'fna_rankC'))
    all_files = [os.path.join(large_data, f) for f in all_files] # Transform to abs path

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
        'csv',
        data_files = file_dict
    )

    dataset = dataset.filter(lambda example: ((example['ranks'] is not None) and (example['__index_level_0__'] is not None)))

    print(dataset)
    print(type(dataset))
    print(dataset['validation']['ranks'][:2])
    print(dataset['validation']['__index_level_0__'][:2])

    # Load in wordlevel tokenizer:
    tokenizer = PreTrainedTokenizerFast(
        model_max_length = 512,
        tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_rankC.json'),
        mask_token = '[MASK]',
        pad_token = '[PAD]',
        cls_token = '[CLS]',
        sep_token = '[SEP]'
    )

    def tokenize_function(examples):
        try:
            return tokenizer(examples['__index_level_0__'], examples['ranks'], truncation=True, padding=True, return_tensors = 'np')
        except:
            print('Problem', examples['__index_level_0__'], examples['ranks'])

    tokenized_dataset = dataset.map(tokenize_function, batched = True, num_proc = 8, remove_columns=['ranks', ' amino', '__index_level_0__'])

    print(tokenized_dataset)
    print(tokenized_dataset['train'][0])

    # Verify that they're the same length:
    i = 0
    lens = []
    for f in tokenized_dataset['train']:
        lens.append(len(f['input_ids']))
        if i > 1000:
            break
        i += 1

    print(sum(lens) / len(lens))

    #exit()

    # def group_texts(examples, block_size = 128):
    #     # Concatenate all texts.
    #     concatenated_examples = {k: sum(examples[k], []) for k in examples.keys()}
    #     total_length = len(concatenated_examples[list(examples.keys())[0]])
    #     # We drop the small remainder, we could add padding if the model supported it instead of this drop, you can
    #         # customize this part to your needs.
    #     total_length = (total_length // block_size) * block_size
    #     # Split by chunks of max_len.
    #     result = {
    #         k: [t[i : i + block_size] for i in range(0, total_length, block_size)]
    #         for k, t in concatenated_examples.items()
    #     }
    #     result["labels"] = result["input_ids"].copy()
    #     return result


    tokenized_dataset.save_to_disk('SavedDatasets/rank_amino.hgf')

    def prep_joint(examples):
        #examples['token_type_ids'] = None
        # Make labels:
        labels = []
        for i in range(len(examples)):
            e = examples['token_type_ids'][i]
            L = [-100] * len(e)
            for j in range(len(e)):
                # Checks if in amino region or the input_id is SEP token:
                if (e[j] == 1) or (examples['input_ids'][j] == 0):
                    break
                L[j] = examples['input_ids'][j]

            labels.append(L)
        
        examples['labels'] = labels
        return examples


    #gtexts = partial(group_texts, block_size = seq_len)

    lm_datasets = tokenized_dataset.map(prep_joint, batched=True, batch_size = batch_size, num_proc=8)

    lm_datasets.save_to_disk('SavedDatasets/rank_amino.hgf')

def from_train():
    tokenized_dataset = load_from_disk('SavedDatasets/rank_amino.hgf')
    
    print(tokenized_dataset)

    def prep_joint(examples):
        #examples['token_type_ids'] = None
       # Make labels:
        labels = []
        e = examples['token_type_ids']
        L = [-100] * len(e)
        for i in range(len(e)):
            #print(e)
            #for j in range(len(e)):
                # Checks if in amino region or the input_id is SEP token:
            if (e[i] == 1) or (examples['input_ids'][i] == 0):
                break
            L[i] = examples['input_ids'][i]

            #labels.append(L)
        
        examples['labels'] = L
        #print('AFTER LABEL ASSIGN')
        return examples


    #gtexts = partial(group_texts, block_size = seq_len)

    lm_datasets = tokenized_dataset.map(prep_joint, batched=False, num_proc = 8)#batch_size = batch_size, num_proc=8)

    lm_datasets.save_to_disk('SavedDatasets/rank_amino.hgf')

if __name__ == '__main__':
    #prep_dataset()
    from_train()
