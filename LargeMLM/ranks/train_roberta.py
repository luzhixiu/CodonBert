import os
import torch
from torch.utils.tensorboard import SummaryWriter

from transformers import Trainer, TrainingArguments, DataCollatorForLanguageModeling
from transformers import RobertaForMaskedLM, RobertaConfig
from transformers import BertTokenizerFast, PreTrainedTokenizerFast
from datasets import load_from_disk

device = "cuda:0" if torch.cuda.is_available() else "cpu"

print('Loading dataset...')
dataset = load_from_disk('SavedDatasets/rank.hgf')
print('Loaded dataset')
dataset.set_format('pt') # Should transition dataset to GPU?
tokenizer = PreTrainedTokenizerFast(
    tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_rank.json')
)

config = RobertaConfig(
    vocab_size = len(tokenizer), # Len gets whole size of vocab (including special tokens)
    torch_dtype = 'float16',
)

model = RobertaForMaskedLM(config).to(device) # Sets model to device

print('Model Set')

print('Num. parameters', model.num_parameters())

collator = DataCollatorForLanguageModeling(tokenizer = tokenizer, mlm = True, mlm_probability=0.15)

training_args = TrainingArguments(
    output_dir = 'Models/model_checkpoints/test',
    overwrite_output_dir = True,
    num_train_epochs = 1, # Set to 1 for now (training)
    learning_rate = 1e-4, 
    per_device_train_batch_size=256, # Following Roberta paper
    save_steps = 1000, # Saves model at every 1000 steps
)

trainer = Trainer(
    model = model,
    args = training_args,
    data_collator = collator,
    train_dataset = dataset['train'],
    eval_dataset = dataset['validation'],
    tb_writer = SummaryWriter(logdir = 'Models/logdirs/test')
)
print('Trainer set\nTraining...')

trainer.train()
trainer.save_model('Models/trained_models/test')