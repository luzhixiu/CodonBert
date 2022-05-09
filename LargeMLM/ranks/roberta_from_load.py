import os
import torch
from torch.utils.tensorboard import SummaryWriter
import math
import time
from tqdm import trange

from transformers import Trainer, TrainingArguments, DataCollatorForLanguageModeling
from transformers import RobertaForMaskedLM, RobertaConfig
from transformers import BertTokenizerFast, PreTrainedTokenizerFast
from transformers.integrations import TensorBoardCallback
from datasets import load_from_disk

device = "cuda:0" if torch.cuda.is_available() else "cpu"
torch.cuda.set_device(device)

print('Loading dataset...')
dataset = load_from_disk('SavedDatasets/rank.hgf')
print('Loaded dataset')
dataset.set_format('pt') # Should transition dataset to GPU?
tokenizer = PreTrainedTokenizerFast(
    tokenizer_file = os.path.join('SavedTokenizers', 'wordpiece_rank.json'),
    mask_token = '[MASK]',
    pad_token = '[PAD]',
    cls_token = '[CLS]',
)#.to(device)

config = RobertaConfig(
    vocab_size = len(tokenizer), # Len gets whole size of vocab (including special tokens)
    max_position_embeddings = 514,
    hidden_size = 256,
    num_hidden_layers = 4,
    num_attention_heads = 4,
    intermediate_size = 1024,
)

model = RobertaForMaskedLM(config)#.to(device) # Sets model to device

model.from_pretrained('Models/trained_models/test').to(device)

print('Model Set')

print('Num. parameters', model.num_parameters())

collator = DataCollatorForLanguageModeling(tokenizer = tokenizer, mlm = True, mlm_probability=0.15)

training_args = TrainingArguments(
    output_dir = 'Models/model_checkpoints/test',
    overwrite_output_dir = True,
    num_train_epochs = 1, # Set to 1 for now (training)
    learning_rate = 1e-4, 
    gradient_accumulation_steps = 1,
    per_device_train_batch_size= 256,
    per_device_eval_batch_size = 8,
    save_steps = 1000, # Saves model at every 1000 steps
    logging_steps = 100,
    fp16=True,
)

trainer = Trainer(
    model = model,
    args = training_args,
    data_collator = collator,
    train_dataset = dataset['train'],
    eval_dataset = dataset['validation'],
    tokenizer = tokenizer,
    callbacks = [TensorBoardCallback(SummaryWriter(log_dir = 'Models/logdirs/test-'+str(int(time.time()))))],
)
print('Trainer set\nTraining...')

# counts = [0] * 10
# size = []

# for i in trange(len(trainer.train_dataset)):
#     for j in trainer.train_dataset[i]['input_ids']:
#         counts[j.item()] += 1
#     size.append(trainer.train_dataset[i]['input_ids'].shape[0])

#     if i > 50_000:
#         break

# print(counts)
# print('Avg size', sum(size) / len(size))

# counts = [0] * 10
# size = []
# for i in trange(len(trainer.eval_dataset)):
#     for j in trainer.eval_dataset[i]['input_ids']:
#         counts[j.item()] += 1
#     size.append(trainer.train_dataset[i]['input_ids'].shape[0])

# print(counts)
# print('Avg size', sum(size) / len(size))
# exit()

#trainer.train()
#trainer.save_model('Models/trained_models/test')


eval_results = trainer.evaluate()
print(f"Loss: {eval_results['eval_loss']};  Perplexity: {math.exp(eval_results['eval_loss']):.2f}")