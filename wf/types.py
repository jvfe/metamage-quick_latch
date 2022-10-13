from dataclasses import dataclass
from enum import Enum

from dataclasses_json import dataclass_json
from latch.types import LatchFile


@dataclass_json
@dataclass
class Sample:
    read1: LatchFile
    read2: LatchFile
    sample_name: str
