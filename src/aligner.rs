use crate::alignment::Alignment;
use crate::config::{AlignmentConfig, AlignmentElement};
use crate::alignment_mtx::{AlignmentMtx, Element};
use crate::alignment_mtx;

pub trait Aligner<S: AlignmentElement, R: AlignmentElement> {
    type Config: AlignmentConfig<S, R>;

    fn align<'a>(&self, subject: &'a str, reference: &'a str, config: &Self::Config) -> Alignment<'a> {
        let mut mtx: AlignmentMtx = self.create_mtx(subject, reference);
        self.fill_top_row(&mut mtx, &config);
        self.fill_left_column(&mut mtx, &config);
        self.fill(&mtx, &config);
        let max: Element = self.find_max(&mtx);
        self.trace_back(&mtx, &max)
    }

    fn create_mtx(&self, subject: &str, reference: &str) -> AlignmentMtx;

    fn fill_top_row(&self, mtx: &mut AlignmentMtx, config: &Self::Config);

    fn fill_left_column(&self, mtx: &mut AlignmentMtx, config: &Self::Config);

    fn fill(&self, mtx: &AlignmentMtx, config: &Self::Config);

    fn find_max(&self, mtx: &AlignmentMtx) -> alignment_mtx::Element;

    fn trace_back<'a>(&self, mtx: &AlignmentMtx, max: &alignment_mtx::Element) -> Alignment<'a>;
}
