//! Reading of variants.

use std::{collections::HashMap, io::BufRead};

use hgvs::static_data::{Assembly, ASSEMBLY_INFOS};
use noodles::vcf::{Header as VcfHeader, Record as VcfRecord};
use noodles_util::variant::reader::{Builder as VariantReaderBuilder, Reader as VariantReader};

/// Provide mapping from contig names to numeric contig IDs.
pub struct ContigMap {
    /// The corresponding assembly.
    pub assembly: Assembly,
    /// Map from contig name to contig ID.
    pub name_map: HashMap<String, usize>,
}

impl ContigMap {
    /// Create a new mapping with the given assembly.
    ///
    /// NB: Grch37 does not include chrMT, Grch7p10 does.
    pub fn new(assembly: Assembly) -> Self {
        let mut name_map = HashMap::new();
        let info = &ASSEMBLY_INFOS[assembly];
        for (idx, seq) in info.sequences.iter().enumerate() {
            name_map.insert(seq.name.clone(), idx);
            for alias in seq.aliases.iter() {
                name_map.insert(alias.clone(), idx);
            }
        }

        Self { assembly, name_map }
    }

    pub fn chrom_to_idx(&self, chrom: &noodles::vcf::record::Chromosome) -> usize {
        match chrom {
            noodles::vcf::record::Chromosome::Name(s)
            | noodles::vcf::record::Chromosome::Symbol(s) => *self
                .name_map
                .get(s)
                .unwrap_or_else(|| panic!("Invalid contig {s}")),
        }
    }
}

/// Read through multiple VCF files at the same time, using
pub struct MultiVcfReader {
    /// The contig mapping to use.
    contig_map: ContigMap,
    /// One reader per file to read from.
    #[allow(clippy::vec_box)]
    readers: Vec<Box<VariantReader<Box<dyn BufRead>>>>,
    /// The headers as read from `readers`.
    #[allow(clippy::vec_box)]
    headers: Vec<Box<VcfHeader>>,
    /// The next record from each reader, if any.
    nexts: Vec<Option<VcfRecord>>,
    /// The smallest from nexts.
    next: usize,
}

impl MultiVcfReader {
    /// Create a new multi VCF reader.
    pub fn new(paths: &[&str], initial_assembly: Option<Assembly>) -> Result<Self, anyhow::Error> {
        let mut assembly: Option<Assembly> = None;

        let mut readers = Vec::new();
        let mut headers = Vec::new();
        let mut nexts = Vec::new();
        for path in paths {
            tracing::trace!("Opening file {}", path);
            let mut reader = Box::new(VariantReaderBuilder::default().build_from_path(path)?);
            let header = Box::new(reader.as_mut().read_header()?);
            assembly = Some(guess_assembly(header.as_ref(), true, initial_assembly)?);
            nexts.push(
                reader
                    .as_mut()
                    .records(header.as_ref())
                    .next()
                    .transpose()?,
            );
            readers.push(reader);
            headers.push(header);
        }

        let contig_map = ContigMap::new(assembly.unwrap());

        Ok(Self {
            contig_map,
            readers,
            headers,
            nexts,
            next: 0,
        })
    }

    /// Obtain the number of input files.
    pub fn input_count(&self) -> usize {
        self.readers.len()
    }

    /// Helper that returns integer pair for `VcfRecord`.
    fn record_to_pair(&self, record: &VcfRecord) -> (usize, usize) {
        let chrom = self.contig_map.chrom_to_idx(record.chromosome());
        let pos: usize = record.position().into();
        (chrom, pos)
    }

    /// Return reference to next record from all input files.
    pub fn peek(&self) -> (&Option<VcfRecord>, usize) {
        (&self.nexts[self.next], self.next)
    }

    /// Pop the next record and read next record from that file into reader.
    pub fn pop(&mut self) -> Result<(Option<VcfRecord>, usize), anyhow::Error> {
        let result = (self.nexts[self.next].clone(), self.next);
        self.nexts[self.next] = self.readers[self.next]
            .as_mut()
            .records(&self.headers[self.next])
            .next()
            .transpose()?;

        self.next = 0;
        for i in 1..self.nexts.len() {
            match (&self.nexts[self.next], &self.nexts[i]) {
                (None, Some(_)) => self.next = i,
                (Some(lhs), Some(rhs)) => {
                    let lhs = self.record_to_pair(lhs);
                    let rhs = self.record_to_pair(rhs);
                    if rhs < lhs {
                        self.next = i;
                    }
                }
                _ => (), // noop
            }
        }

        Ok(result)
    }
}

/// Canonical chromosome names.
///
/// Note that the mitochondrial genome runs under two names.
const CANONICAL: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "M", "MT",
];

/// Guess the assembly from the given header.
///
/// If the header only contains chrM, for example, the result may be ambiguous. Use `ambiguous_ok`
/// to allow or disallow this.  You can specify an initial value for the assembly to overcome
/// issues.  If the result is incompatible with the `initial_assembly` then an error will
/// be returned.
pub fn guess_assembly(
    vcf_header: &VcfHeader,
    ambiguous_ok: bool,
    initial_assembly: Option<Assembly>,
) -> Result<Assembly, anyhow::Error> {
    let mut result = initial_assembly;

    for (assembly, info) in ASSEMBLY_INFOS.iter() {
        let contig_map = ContigMap::new(assembly);
        let mut lengths = HashMap::new();
        for seq in &info.sequences {
            if CANONICAL.contains(&seq.name.as_str()) {
                lengths.insert(
                    contig_map.name_map.get(seq.name.as_str()).unwrap(),
                    seq.length,
                );
            }
        }

        let mut incompatible = 0;
        let mut compatible = 0;
        for (name, data) in vcf_header.contigs() {
            if let Some(length) = data.length() {
                let idx = contig_map.name_map.get(name.as_ref());
                if let Some(idx) = idx {
                    let name = &info.sequences[*idx].name;
                    if CANONICAL.contains(&name.as_ref()) {
                        if *lengths.get(idx).unwrap() == length {
                            compatible += 1;
                        } else {
                            incompatible += 1;
                        }
                    }
                }
            } else {
                tracing::warn!(
                    "Cannot guess assembly because no length for contig {}",
                    &name
                );
                compatible = 0;
                break;
            }
        }

        if compatible > 0 && incompatible == 0 {
            if let Some(value) = result {
                if !ambiguous_ok {
                    return Err(anyhow::anyhow!(
                        "Found two matching assemblies {:?} / {:?}",
                        value,
                        assembly
                    ));
                } else if result != initial_assembly {
                    result = Some(assembly);
                }
            } else {
                result = Some(assembly);
            }
        } else if initial_assembly.is_some() && assembly == initial_assembly.unwrap() {
            return Err(anyhow::anyhow!(
                "Incompatible with initial assembly {:?}",
                result.unwrap()
            ));
        }
    }

    if let Some(result) = result {
        Ok(result)
    } else {
        Err(anyhow::anyhow!("No matching assembly found"))
    }
}

#[cfg(test)]
mod test {
    use hgvs::static_data::Assembly;

    use super::*;
    use noodles_util::variant::reader::Builder as VariantReaderBuilder;

    #[test]
    fn guess_assembly_helix_chrmt_ambiguous_ok_initial_none() -> Result<(), anyhow::Error> {
        let path = "tests/data/db/create/seqvar_freqs/helix.chrM.vcf";
        let mut reader = VariantReaderBuilder::default().build_from_path(path)?;
        let header = reader.read_header()?;

        let actual = guess_assembly(&header, true, None)?;
        assert_eq!(actual, Assembly::Grch38);

        Ok(())
    }

    #[test]
    fn guess_assembly_helix_chrmt_ambiguous_ok_initial_override() -> Result<(), anyhow::Error> {
        let path = "tests/data/db/create/seqvar_freqs/helix.chrM.vcf";
        let mut reader = VariantReaderBuilder::default().build_from_path(path)?;
        let header = reader.read_header()?;

        let actual = guess_assembly(&header, true, Some(Assembly::Grch37p10))?;
        assert_eq!(actual, Assembly::Grch37p10);

        Ok(())
    }

    #[test]
    fn guess_assembly_helix_chrmt_ambiguous_ok_initial_override_fails() -> Result<(), anyhow::Error>
    {
        let path = "tests/data/db/create/seqvar_freqs/helix.chrM.vcf";
        let mut reader = VariantReaderBuilder::default().build_from_path(path)?;
        let header = reader.read_header()?;

        assert!(guess_assembly(&header, true, Some(Assembly::Grch37)).is_err());

        Ok(())
    }

    #[test]
    fn guess_assembly_helix_chrmt_ambiguous_fail() -> Result<(), anyhow::Error> {
        let path = "tests/data/db/create/seqvar_freqs/helix.chrM.vcf";
        let mut reader = VariantReaderBuilder::default().build_from_path(path)?;
        let header = reader.read_header()?;

        assert!(guess_assembly(&header, false, None).is_err());

        Ok(())
    }

    #[test]
    fn contig_map_smoke() {
        ContigMap::new(Assembly::Grch37p10);
        ContigMap::new(Assembly::Grch38);
    }

    #[test]
    fn test_multivcf_reader() -> Result<(), anyhow::Error> {
        let mut reader = MultiVcfReader::new(
            &[
                "tests/data/db/create/seqvar_freqs/mt/gnomad.chrM.vcf",
                "tests/data/db/create/seqvar_freqs/mt/helix.chrM.vcf",
            ],
            Some(Assembly::Grch37p10),
        )?;

        let top = reader.peek();
        assert_eq!(top.1, 0);
        if let Some(record) = top.0 {
            assert_eq!(Into::<usize>::into(record.position()), 3);
        };

        for _i in 0..6 {
            assert!(reader.pop()?.0.is_some());
        }
        assert!(reader.peek().0.is_none());

        Ok(())
    }
}
