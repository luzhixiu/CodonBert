import os, random
import torch
import pandas as pd
from datasets import load_dataset
from transformers import Trainer, TrainingArguments
from transformers import BigBirdModel
from transformers import BertTokenizerFast, PreTrainedTokenizerFast

# Change this for other users:
owens_desktop = '/Users/owenqueen/Desktop/bioinformatics/codonbert/CodonBert/Data'
prefix = owens_desktop

def prep_dataset():
    large_data = os.path.join(prefix, 'fnaCollection')
    all_files = os.listdir(os.path.join(prefix, 'fnaCollection'))
    all_files = [os.path.join(large_data, f) for f in all_files] # Transform to abs path

    val_pct = 0.1

    # Choose val_pct of all_files for validation data:
    val_sample = random.sample(range(len(all_files)), k = int(len(all_files) * val_pct))
    val_files = [all_files[i] for i in val_sample]
    val_set = set(val_sample) # Convert to set for O(1) lookup
    train_files = [all_files[i] for i in range(len(all_files)) if i not in val_set]

    file_dict = {
        'train': train_files,
        'validation': val_files,
        'test': os.path.join(prefix, 's288c.fasta.test.txt') # Test on predefined testing file
    }

    dataset = load_dataset(
        'text',
        data_files = file_dict
    )

    print(dataset)
    print(type(dataset))

    # Load in wordlevel tokenizer:
    tokenizer = PreTrainedTokenizerFast(tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_codon.json'))

    def tokenize_function(examples):
        return tokenizer(examples['text'])

    tokenized_dataset = dataset.map(tokenize_function, batched = True, num_proc = 4, remove_columns=['text'])

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

    lm_datasets = tokenized_dataset.map(group_texts, batched=True, batch_size = 1000, num_proc=4)

    # print(tokenized_datasets)
    # print(tokenized_datasets['train'][0])

if __name__ == '__main__':
    prep_dataset()
