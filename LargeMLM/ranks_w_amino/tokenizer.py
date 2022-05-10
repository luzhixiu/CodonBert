import os
from tokenizers import decoders, models, normalizers, pre_tokenizers, processors, trainers, Tokenizer
from datasets import load_dataset
    
prefix = '/lustre/isaac/scratch/oqueen/CodonBert/Data'
large_data = os.path.join(prefix, 'fna_rankC')
all_files = os.listdir(os.path.join(prefix, 'fna_rankC'))
all_files = [os.path.join(large_data, f) for f in all_files]

codonTable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def train_wordpiece_tokenizer():

    # Make dataset:
    owens_desktop = '/Users/owenqueen/Desktop/bioinformatics/codonbert/CodonBert/Data'
    prefix = owens_desktop

    # splits = []
    # for f in all_files:
    #     L = open(f, 'r').readlines()
    #     splits += L

    #splits = [['']]

    #print(len(splits))

    # Make trivial length of strings with all numbers:
    batch_iterator = [str(i) + ' ' + str(i+1) for i in range(1, 6)]
    batch_iterator += [f'{k} {k}' for k in codonTable.values()]
    #print(batch_iterator)

    tokenizer = gen_tokenizer(batch_iterator=batch_iterator)

    tokenizer.save(os.path.join('SavedTokenizers', 'wordpiece_rankC.json'))


def gen_tokenizer(vocab_json_path = None, batch_iterator = None):
    tokenizer = Tokenizer(models.WordLevel(unk_token='[UNK]'))

    if vocab_json_path is None:
        tokenizer.normalizer = normalizers.BertNormalizer(lowercase=False)

        tokenizer.normalizer = normalizers.Sequence(
            [normalizers.NFD(),normalizers.StripAccents()]
        )
        tokenizer.pre_tokenizer = pre_tokenizers.BertPreTokenizer()

        # Will eventually add species-level special tokens:
        special_tokens = ["[UNK]", "[PAD]", "[CLS]", "[MASK]"]
        #trainer = trainers.WordPieceTrainer(vocab_size=5000, special_tokens=special_tokens)
        trainer = trainers.WordLevelTrainer(special_tokens = special_tokens)
        tokenizer.train_from_iterator(batch_iterator, trainer = trainer)

        # Post processing for BERT: https://huggingface.co/docs/tokenizers/python/latest/pipeline.html
        tokenizer.post_processor = processors.TemplateProcessing(
            single="[CLS] $0 [SEP]",
            pair="[CLS] $A [SEP] $B:1 [SEP]:1",
            special_tokens=[("[CLS]", 1), ("[SEP]", 0)],
        )

    else:
        # Load from vocab.json
        pass

    return tokenizer

if __name__ == '__main__':
    train_wordpiece_tokenizer()