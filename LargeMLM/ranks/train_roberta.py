import os
import torch

from transformers import Trainer, TrainingArguments, DataCollatorForLanguageModeling
from transformers import RobertaForMaskedLM, RobertaConfig
from transformers import BertTokenizerFast, PreTrainedTokenizerFast
from datasets import load_from_disk

device = "cuda:0" if torch.cuda.is_available() else "cpu"

dataset = load_from_disk('SavedDatasets/rank.hgf')
tokenizer = PreTrainedTokenizerFast(tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_rank.json'))

config = RobertaConfig(
    vocab_size = len(tokenizer)
)

model = RobertaForMaskedLM(config).to(device)

print('Num. parameters', model.num_parameters())

collator = DataCollatorForLanguageModeling(tokenizer = tokenizer, mlm = True, mlm_probability=0.15)

training_args = TrainingArguments(
    output_dir = 'trained_models/test',
    overwrite_output_dir = True,
    num_train_epochs = 1,
    learning_rate = 1e-4,
    per_device_train_batch_size=256,
    save_steps = 1000,
)