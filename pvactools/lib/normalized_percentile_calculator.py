import os
import threading
import requests
import h5py
import numpy as np

class PercentileFileCache:
    """Thread-safe cache for allele file download status."""
    _missing = set()
    _events = {}
    _lock = threading.Lock()

    @classmethod
    def get_event(cls, allele):
        with cls._lock:
            if allele not in cls._events:
                cls._events[allele] = threading.Event()
            return cls._events[allele]

    @classmethod
    def mark_missing(cls, allele):
        with cls._lock:
            cls._missing.add(allele)

    @classmethod
    def is_missing(cls, allele):
        with cls._lock:
            return allele in cls._missing

    @classmethod
    def mark_done(cls, allele):
        event = cls.get_event(allele)
        event.set()  # unblock all waiting threads

class NormalizedPercentileCalculator:
    def __init__(self, reference_scores_path):
        self.reference_scores = {}
        self.reference_scores_path = reference_scores_path

    def calculate_normalized_percentile(self, allele, length, score, method, is_reversed=False, mode="length_agnostic"):
        if allele is None or length is None or score is None or score == 'NA':
            return 'NA'

        normalized = self.normalize_allele(allele)
        if normalized is None:
            return 'NA'

        if not self.has_normalized_file(allele):
            return 'NA'

        allele_file = f"{normalized}_percentiles.h5"
        file_path = os.path.join(self.reference_scores_path, allele_file)

        if mode == "per_length":
            key = f"{method}/{length}mer"
        elif mode == "length_agnostic":
            key = method
        cache_key = f"{normalized}_{key}"
        if cache_key in self.reference_scores:
            ref_scores = self.reference_scores[cache_key]
        else:
            try:
                with h5py.File(file_path, "r") as f:
                    if key not in f:
                        return 'NA'  # algorithm or length not present
                    if mode == "per_length":
                        ref_scores = f[key][...]
                    elif mode == "length_agnostic":
                        scores = []
                        for length in f[key].keys():
                            scores.extend(f[key][length][...])
                        ref_scores = np.array(sorted(scores))
                    self.reference_scores[cache_key] = ref_scores
            except Exception as e:
                return 'NA'

        if ref_scores.size == 0:
            return 'NA'

        n = len(ref_scores)
        left = np.searchsorted(ref_scores, score, side="left")
        right = np.searchsorted(ref_scores, score, side="right")

        if left == right:
            percentile = left / n * 100
        else:
            percentile = (left + right) / (2 * n) * 100

        percentile = 100 - percentile if is_reversed else percentile
        return round(percentile, 3)

    def has_normalized_file(self, allele):
        normalized_allele = self.normalize_allele(allele)
        allele_file = f"{normalized_allele}_percentiles.h5"
        file_path = os.path.join(self.reference_scores_path, allele_file)

        if os.path.exists(file_path):
            return True

        event = PercentileFileCache.get_event(normalized_allele)

        if PercentileFileCache.is_missing(normalized_allele) and event.is_set():
            return False

        with PercentileFileCache._lock:
            if os.path.exists(file_path):
                return True

            if not event.is_set() and normalized_allele not in PercentileFileCache._missing:
                pass
            else:
                # Another thread is already downloading – wait for it
                event.wait()
                return os.path.exists(file_path)

        url = f"https://raw.githubusercontent.com/griffithlab/pvactools_percentiles_data/main/hdf5/{allele_file}"
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()

            os.makedirs(self.reference_scores_path, exist_ok=True)
            with open(file_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            PercentileFileCache.mark_done(normalized_allele)
            return True

        except Exception:
            PercentileFileCache.mark_missing(normalized_allele)
            PercentileFileCache.mark_done(normalized_allele)

            print(f"WARNING: No percentile file found for allele {allele} ({normalized_allele}). Using 'NA' for percentile values.")
            return False

    def normalize_allele(self, allele):
        if allele is None:
            return None

        return allele.replace("*", "_").replace(":", "_")
