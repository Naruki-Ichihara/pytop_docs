import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Extended FEniCS for Topology Optimization',
    Svg: require('@site/static/img/AcademicCap.svg').default,
    description: (
      <>
        Pytop extends the capabilities of FEniCS specifically for topology optimization, providing robust tools for advanced computational design.
      </>
    ),
  },
  {
    title: 'MPI Ready',
    Svg: require('@site/static/img/ServerStack.svg').default,
    description: (
      <>
        Pytop is fully compatible with MPI, allowing you to monitor parallel applications efficiently.
      </>
    ),
  },
  {
    title: 'Multiphysics Support',
    Svg: require('@site/static/img/Bolt.svg').default,
    description: (
      <>
        Pytop supports multiphysics simulations, integrating seamlessly with FEniCS to handle complex interactions between different physical phenomena.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} role="img" />
      </div>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
