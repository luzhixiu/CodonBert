import os
from tokenizers import decoders, models, normalizers, pre_tokenizers, processors, trainers, Tokenizer
from datasets import load_dataset
    

def train_wordpiece_tokenizer():

    # Make dataset:
    owens_desktop = '/Users/owenqueen/Desktop/bioinformatics/codonbert/CodonBert/Data'
    prefix = owens_desktop

    # Make trivial length of strings with all numbers:
    batch_iterator = [str(i) + ' ' + str(i+1) for i in range(1, 6)]
    print(batch_iterator)

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
        special_tokens = ["[UNK]", "[PAD]", "[CLS]", "[MASK]"]
        #trainer = trainers.WordPieceTrainer(vocab_size=5000, special_tokens=special_tokens)
        trainer = trainers.WordLevelTrainer(special_tokens = special_tokens)
        tokenizer.train_from_iterator(batch_iterator, trainer = trainer)

        # Post processing for BERT: https://huggingface.co/docs/tokenizers/python/latest/pipeline.html
        tokenizer.post_processor = processors.TemplateProcessing(
            single="[CLS] $A",
            special_tokens=[("[CLS]", 1)],
        )

    else:
        # Load from vocab.json
        pass

    return tokenizer

if __name__ == '__main__':
    train_wordpiece_tokenizer()