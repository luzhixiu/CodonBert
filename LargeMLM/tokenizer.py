import os
from tokenizers import decoders, models, normalizers, pre_tokenizers, processors, trainers, Tokenizer
from datasets import load_dataset


def make_codon_int_mapping():
    codon_int_mapping = {
    'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
    'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,'ATA': 0, 'ATG': 0,
    'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
    'TAA': 0, 'TAG': 0,'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
    'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
    'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,   'TCG': 0,
    'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,'ACT': 0, 'ACC': 0,
    'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
    'TGT': 0, 'TGC': 0,'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
    'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
    }

    i = 0
    for k in codon_int_mapping.keys():
        codon_int_mapping[k] = i
        i += 1
    

def train_wordpiece_tokenizer():
    # codon_int_mapping = {
    #     'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0,
    #     'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,'ATA': 0, 'ATG': 0,
    #     'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'TAT': 0, 'TAC': 0,
    #     'TAA': 0, 'TAG': 0,'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0,
    #     'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
    #     'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,   'TCG': 0,
    #     'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,'ACT': 0, 'ACC': 0,
    #     'ACA': 0, 'ACG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
    #     'TGT': 0, 'TGC': 0,'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
    #     'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
    #     'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
    # }

    # Make dataset:
    owens_desktop = '/Users/owenqueen/Desktop/bioinformatics/codonbert/CodonBert/Data'
    prefix = owens_desktop

    dataset = load_dataset(
        'text',
        data_files = os.path.join(prefix, 'random.train.txt'),
        split='train'
    )
    batch_size = 1000
    all_texts = [dataset[i : i + batch_size]["text"] for i in range(0, len(dataset), batch_size)]

    def batch_iterator():
        for i in range(0, len(dataset), batch_size):
            yield dataset[i : i + batch_size]["text"]

    tokenizer = gen_tokenizer(batch_iterator=batch_iterator)

    tokenizer.save(os.path.join('SavedTokenizers', 'wordpiece_codon.json'))


def gen_tokenizer(vocab_json_path = None, batch_iterator = None):
    tokenizer = Tokenizer(models.WordLevel(unk_token='[UNK]'))

    if vocab_json_path is None:
        tokenizer.normalizer = normalizers.BertNormalizer(lowercase=False)

        tokenizer.normalizer = normalizers.Sequence(
            [normalizers.NFD(),normalizers.StripAccents()]
        )
        tokenizer.pre_tokenizer = pre_tokenizers.BertPreTokenizer()

        # Will eventually add species-level special tokens:
        special_tokens = ["[UNK]", "[PAD]", "[CLS]", "[SEP]", "[MASK]"]
        #trainer = trainers.WordPieceTrainer(vocab_size=5000, special_tokens=special_tokens)
        trainer = trainers.WordLevelTrainer(special_tokens = special_tokens)
        tokenizer.train_from_iterator(batch_iterator(), trainer = trainer)

        # Post processing for BERT: https://huggingface.co/docs/tokenizers/python/latest/pipeline.html
        tokenizer.post_processor = processors.TemplateProcessing(
            single="[CLS] $A [SEP]",
            pair="[CLS] $A [SEP] $B:1 [SEP]:1",
            special_tokens=[("[CLS]", 1), ("[SEP]", 2)],
        )

    else:
        # Load from vocab.json
        pass

    return tokenizer

if __name__ == '__main__':
    train_wordpiece_tokenizer()