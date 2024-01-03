#[derive(Clone)]
pub enum Duration {
    Days(f64),
    Seconds(f64),
    Years(f64),
    Minutes(f64),
    Hours(f64),
}

impl Duration {
    pub fn days(&self) -> f64 {
        match self {
            Duration::Days(v) => *v,
            Duration::Seconds(v) => *v / 86400.0,
            Duration::Years(v) => *v * 365.25,
            Duration::Minutes(v) => *v / 1440.0,
            Duration::Hours(v) => *v / 24.0,
        }
    }
    pub fn seconds(&self) -> f64 {
        match self {
            Duration::Days(v) => *v * 86400.0,
            Duration::Seconds(v) => *v,
            Duration::Years(v) => *v * 86400.0 * 365.25,
            Duration::Minutes(v) => *v * 60.0,
            Duration::Hours(v) => *v * 3600.0,
        }
    }
    pub fn hours(&self) -> f64 {
        match self {
            Duration::Days(v) => *v * 24.0,
            Duration::Seconds(v) => *v / 3600.0,
            Duration::Minutes(v) => *v / 60.0,
            Duration::Hours(v) => *v,
            Duration::Years(v) => *v * 24.0 * 365.25,
        }
    }

    pub fn minutes(&self) -> f64 {
        match self {
            Duration::Days(v) => *v * 1440.0,
            Duration::Seconds(v) => *v / 60.0,
            Duration::Minutes(v) => *v,
            Duration::Hours(v) => *v * 60.0,
            Duration::Years(v) => *v * 1440.0 * 365.25,
        }
    }

    pub fn to_string(&self) -> String {
        let mut secs = self.seconds();

        if secs < 1.0 {
            format!("Duration: {:.3} microseconds", (secs % 1.0) * 1.0e6)
        } else {
            let days = (secs / 86400.0) as usize;
            let hours = ((secs % 86400.0) / 3600.0) as usize;
            let minutes = ((secs % 3600.0) / 60.0) as usize;
            secs = secs % 60.0;
            let mut s = String::from("Duration: ");
            if days > 0 {
                s.push_str(format!("{} days, ", days).as_str());
            }
            if hours > 0 || days > 0 {
                s.push_str(format!("{} hours, ", hours).as_str());
            }
            if minutes > 0 || hours > 0 || days > 0 {
                s.push_str(format!("{} minutes, ", minutes).as_str());
            }
            s.push_str(format!("{:.3} seconds", secs).as_str());
            s
        }
    }
}

impl std::ops::Add<Duration> for Duration {
    type Output = Duration;
    fn add(self, other: Duration) -> Self::Output {
        Duration::Seconds(self.seconds() + other.seconds())
    }
}

impl std::ops::Add<&Duration> for Duration {
    type Output = Duration;
    fn add(self, other: &Duration) -> Self::Output {
        Duration::Seconds(self.seconds() + other.seconds())
    }
}

impl std::ops::Add<Duration> for &Duration {
    type Output = Duration;
    fn add(self, other: Duration) -> Self::Output {
        Duration::Seconds(self.seconds() + other.seconds())
    }
}

impl std::ops::Sub<Duration> for Duration {
    type Output = Duration;
    fn sub(self, other: Duration) -> Self::Output {
        Duration::Seconds(self.seconds() - other.seconds())
    }
}

impl std::ops::Sub<&Duration> for Duration {
    type Output = Duration;
    fn sub(self, other: &Duration) -> Self::Output {
        Duration::Seconds(self.seconds() - other.seconds())
    }
}

impl std::ops::Sub<Duration> for &Duration {
    type Output = Duration;
    fn sub(self, other: Duration) -> Self::Output {
        Duration::Seconds(self.seconds() - other.seconds())
    }
}

impl std::ops::Sub<&Duration> for &Duration {
    type Output = Duration;
    fn sub(self, other: &Duration) -> Self::Output {
        Duration::Seconds(self.seconds() - other.seconds())
    }
}

impl std::fmt::Display for Duration {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}
