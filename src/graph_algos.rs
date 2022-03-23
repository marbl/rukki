pub mod scc;
pub mod superbubble;
pub mod dfs;

pub fn only_or_none<T>(mut iter: impl Iterator<Item=T>) -> Option<T> {
    let e = iter.next()?;
    match iter.next() {
        None => Some(e),
        _ => None,
    }
}