use num::cast::AsPrimitive;
use std::io::{Write, Cursor};


pub mod filter;
pub mod convert;
pub mod config;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}




#[derive(Copy, Clone)]
pub enum RADType {
    BOOL,
    U8,
    U16,
    U32,
    U64,
    F32,
    F64,
}



pub fn encode_type_tag(type_tag: RADType) -> Option<u8> {
    match type_tag {
        RADType::BOOL => Some(0),
        RADType::U8 => Some(1),
        RADType::U16 => Some(2),
        RADType::U32 => Some(3),
        RADType::U64 => Some(4),
        RADType::F32 => Some(5),
        RADType::F64 => Some(6),
        //_ => None,
    }
}

#[derive(Copy, Clone)]
pub enum RADIntID {
    U8,
    U16,
    U32,
    U64,
}

impl RADIntID {
    pub fn bytes_for_type(&self) -> usize {
        match self {
            Self::U8 => std::mem::size_of::<u8>(),
            Self::U16 => std::mem::size_of::<u16>(),
            Self::U32 => std::mem::size_of::<u32>(),
            Self::U64 => std::mem::size_of::<u64>(),
        }
    }

    /// Based on the variant of the current enum, write the value `v`  
    /// out using `owrite`.  Here, `v` is bound to be some primitive
    /// integer type.  It is the responsibility of the caller to ensure
    /// that, if `v` is wider than the enum type on which this function
    /// is called, no important information is lost by discarding the higher
    /// order bits.
    pub fn write_to<T: PrimitiveInteger, U: Write>(
        &self,
        v: T,
        owriter: &mut U,
    ) -> std::io::Result<()> {
        match self {
            Self::U8 => {
                let vo: u8 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
            Self::U16 => {
                let vo: u16 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
            Self::U32 => {
                let vo: u32 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
            Self::U64 => {
                let vo: u64 = v.as_();
                owriter.write_all(&vo.to_le_bytes())
            }
        }
    }
}




pub trait PrimitiveInteger:
    AsPrimitive<u8>
    + AsPrimitive<u16>
    + AsPrimitive<u32>
    + AsPrimitive<u64>
    + AsPrimitive<usize>
    + AsPrimitive<i8>
    + AsPrimitive<i16>
    + AsPrimitive<i32>
    + AsPrimitive<i64>
    + AsPrimitive<isize>
{
}

impl<
        T: AsPrimitive<u8>
            + AsPrimitive<u16>
            + AsPrimitive<u32>
            + AsPrimitive<u64>
            + AsPrimitive<usize>
            + AsPrimitive<i8>
            + AsPrimitive<i16>
            + AsPrimitive<i32>
            + AsPrimitive<i64>
            + AsPrimitive<isize>,
    > PrimitiveInteger for T
{
}


pub fn write_str_bin(v: &str, type_id: &RADIntID, owriter: &mut Cursor<Vec<u8>>) {
    match type_id {
        RADIntID::U8 => {
            owriter
                .write_all(&(v.len() as u8).to_le_bytes())
                .expect("coudn't write to output file");
        }
        RADIntID::U16 => {
            owriter
                .write_all(&(v.len() as u16).to_le_bytes())
                .expect("coudn't write to output file");
        }
        RADIntID::U32 => {
            owriter
                .write_all(&(v.len() as u32).to_le_bytes())
                .expect("coudn't write to output file");
        }
        RADIntID::U64 => {
            owriter
                .write_all(&(v.len() as u64).to_le_bytes())
                .expect("coudn't write to output file");
        }
    }
    owriter
        .write_all(v.as_bytes())
        .expect("coudn't write to output file");
}