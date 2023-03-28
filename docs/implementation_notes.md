Implementation notes.

## Frequency Databases

Frequency databases are implemented using [RocksDB](https://rocksdb.org/).
There are four column families:

- `meta` - store meta information
- `autosomal` - store autosomal variant frequencies
- `gonosomal` - store gonomosal variant frequencies
- `mitochondrial` - store mitochondrial variant frequencies

The variant is encoded as the key.
See `impl From<Var> for Vec<u8>` on how the serialization takes place.

Currently, we create the RocksDB database without any particular tuning to performance or compression settings.

See the data structures in [`crate::db::create::seqvar_freqs::serialized`] for the serialization of counts as little endian 32 bit integers.

## Transcript Databases

* Transcript databases are stored as [Flatbuffers](https://github.com/google/flatbuffers).
* Array-backed interval trees from [rust-bio](https://github.com/rust-bio/rust-bio) are used for fast lookup from chromosomal coordinate to transcript.
* Transcripts are taken from [cdot](https://github.com/SACGF/cdot).
