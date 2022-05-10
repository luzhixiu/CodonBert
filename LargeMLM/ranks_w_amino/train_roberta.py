import os
#import torch
from torch.utils.tensorboard import SummaryWriter
import math
import time

from transformers import Trainer, TrainingArguments, DataCollatorForLanguageModeling
from transformers import RobertaForMaskedLM, RobertaConfig
from transformers import BertTokenizerFast, PreTrainedTokenizerFast
from transformers.integrations import TensorBoardCallback
from datasets import load_from_disk

#device = "cuda:0" if torch.cuda.is_available() else "cpu"

print('Loading dataset...')
dataset = load_from_disk('SavedDatasets/rank_amino.hgf')
print('Loaded dataset')
#dataset.set_format('pt') # Should transition dataset to GPU?
tokenizer = PreTrainedTokenizerFast(
    model_max_length = 512,
    tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_rankC.json'),
    mask_token = '[MASK]',
    pad_token = '[PAD]',
    cls_token = '[CLS]',
    sep_token = '[SEP]'
)

config = RobertaConfig(
    vocab_size = len(tokenizer), # Len gets whole size of vocab (including special tokens)
    max_position_embeddings = 1000,
    hidden_size = 256,
    num_hidden_layers = 4,
    num_attention_heads = 4,
    intermediate_size = 1024,
)

model = RobertaForMaskedLM(config)#.to(device) # Sets model to device

print('Model Set')

print('Num. parameters', model.num_parameters())

collator = DataCollatorForLanguageModeling(tokenizer = tokenizer, mlm = True, mlm_probability=0.15)

training_args = TrainingArguments(
    output_dir = 'Models/model_checkpoints/rankC_test',
    overwrite_output_dir = True,
    num_train_epochs = 10, # Set to 1 for now (training)
    learning_rate = 1e-4, 
    gradient_accumulation_steps = 1,
    per_device_train_batch_size= 256,
    per_device_eval_batch_size = 8,
    save_steps = 1000, # Saves model at every 1000 steps
    logging_steps = 100,
    evaluation_strategy = 'no',
    fp16=True,
)

trainer = Trainer(
    model = model,
    args = training_args,
    data_collator = collator,
    train_dataset = dataset['train'],
    eval_dataset = dataset['validation'],
    tokenizer = tokenizer,
    callbacks = [TensorBoardCallback(SummaryWriter(log_dir = 'Models/logdirs/rankC-test-'+str(int(time.time()))))],
)
print('Trainer set\nTraining...')

trainer.train()
trainer.save_model('Models/trained_models/rankC_test')


eval_results = trainer.evaluate()
print(f"Loss: {eval_results['eval_loss']};  Perplexity: {math.exp(eval_results['eval_loss']):.2f}")


