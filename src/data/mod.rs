//! Data preparation and aggregation.

pub mod aggregate;
pub mod soa;

pub use aggregate::{AggregatedMatrix, AggregationMode, ClusterMap, aggregate, load_cluster_map};
pub use soa::{
    ExpressionSoA, PreparedExpression, prepare_expression, prepare_expression_with_clusters,
};
