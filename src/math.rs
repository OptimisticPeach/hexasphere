#[cfg(feature = "libm")]
mod libm_math {
    #[inline(always)]
    pub(crate) fn acos(f: f32) -> f32 {
        libm::acosf(f)
    }

    #[allow(unused)]
    #[inline(always)]
    pub(crate) fn sin(f: f32) -> f32 {
        libm::sinf(f)
    }

    #[inline(always)]
    pub(crate) fn sqrt(f: f32) -> f32 {
        libm::sqrtf(f)
    }
}

#[cfg(not(feature = "libm"))]
mod std_math {
    #[inline(always)]
    pub(crate) fn acos(f: f32) -> f32 {
        f32::acos(f)
    }

    #[allow(unused)]
    #[inline(always)]
    pub(crate) fn sin(f: f32) -> f32 {
        f32::sin(f)
    }

    #[inline(always)]
    pub(crate) fn sqrt(f: f32) -> f32 {
        f32::sqrt(f)
    }
}

#[cfg(feature = "libm")]
pub(crate) use libm_math::*;

#[cfg(not(feature = "libm"))]
pub(crate) use std_math::*;
