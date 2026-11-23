from dataclasses import dataclass
from pathlib import Path


class classproperty(property):
    def __get__(self, cls, owner):
        return classmethod(self.fget).__get__(None, owner)()


@dataclass
class GenomadData:
    data_dir: Path = Path(__file__).parents[0] / "data"

    @classproperty
    def host_nn_model_file(cls) -> Path:
        return Path("/data/DeepMetaClass/training_checkpoints/classifier_undersampling_nb_stacks_2") / "final_classifier_weights.h5"


@dataclass
class GenomadOutputs:
    prefix: str
    output_dir: Path
    
    # host-classification

    @property
    def host_classification_log(self) -> Path:
        return self.output_dir / f"{self.prefix}_host_classification.log"

    @property
    def host_classification_dir(self) -> Path:
        return self.output_dir / f"{self.prefix}_host_classification"

    @property
    def host_classification_execution_info(self) -> Path:
        return self.host_classification_dir / f"{self.prefix}_host_classification.json"

    @property
    def host_encoded_sequences_dir(self) -> Path:
        return self.host_classification_dir / f"{self.prefix}_encoded_sequences"

    @property
    def host_seq_window_id_output(self) -> Path:
        return self.host_encoded_sequences_dir / f"{self.prefix}_seq_window_id.npz"

    @property
    def host_classification_output(self) -> Path:
        return self.host_classification_dir / f"{self.prefix}_host_classification.tsv"

    @property
    def host_classification_npz_output(self) -> Path:
        return self.host_classification_dir / f"{self.prefix}_host_classification.npz"

    @property
    def host_encoded_proviruses_dir(self) -> Path:
        return self.host_classification_dir / f"{self.prefix}_encoded_proviruses"

    @property
    def host_provirus_window_id_output(self) -> Path:
        return self.host_encoded_proviruses_dir / f"{self.prefix}_provirus_window_id.npz"

    @property
    def provirus_host_classification_output(self) -> Path:
        return (self.host_classification_dir / f"{self.prefix}_provirus_host_classification.tsv")

    @property
    def provirus_host_classification_npz_output(self) -> Path:
        return (self.host_classification_dir / f"{self.prefix}_provirus_host_classification.npz")
