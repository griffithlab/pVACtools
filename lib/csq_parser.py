import re

class CsqParser:
    def __init__(self, csq_header_description):
        format_pattern = re.compile('Format: (.*)')
        match = format_pattern.search(csq_header_description)
        self.csq_format = match.group(1)

    def parse_csq_entries_for_allele(self, csq_entries, csq_allele):
        csq_format_array = self.csq_format.split('|')

        transcripts = []
        for entry in csq_entries:
            values = entry.split('|')
            transcript = {}
            for key, value in zip(csq_format_array, values):
                transcript[key] = value
            if transcript['Allele'] == csq_allele:
                transcripts.append(transcript)

        return transcripts

    def resolve_alleles(self, entry):
        alleles = {}
        for alt in entry.ALT:
            alt = str(alt)
            if self.is_insertion(entry.REF, alt) or self.is_deletion(entry.REF, alt):
                if alt[0:1] != entry.REF[0:1]:
                    alleles[alt] = alt
                elif alt[1:] == "":
                    alleles[alt] = '-'
                else:
                    alleles[alt] = alt[1:]
            else:
                alt = str(alt)
                alleles[alt] = alt

        return alleles

    def is_insertion(self, ref, alt):
        return len(alt) > len(ref)

    def is_deletion(self, ref, alt):
        return len(alt) < len(ref)
